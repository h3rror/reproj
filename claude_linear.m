% advdiff_opinf_experiment.m
%
% Numerical check of the exact-recovery theorem for operator inference (OpInf).
%
% Full model: 1D advection-diffusion  u_t + c u_x = nu u_xx  on the periodic unit
% interval, central finite differences in space, FORWARD EULER in time:
%
%     x_{k+1} = M x_k,   M = I + dt*A.
%
% Initial condition: square pulse (discontinuous, i.e. non-smooth).
%
% Theorem: if the POD basis V is computed from the same noise-free snapshot matrix
% X = [x_0,...,x_{K-1}] that is used as regression input, then unregularized OpInf
%     min_Mt || Mt*U - Up ||_F,   U = V'*X,  Up = V'*X_+
% returns exactly  Mhat = V'*M*V  (the intrinsic Galerkin operator).
%
% Control cases break one assumption each and should NOT recover V'*M*V.
%
% Runs in MATLAB and in GNU Octave (no local functions are used).

clear; close all; clc;
EPS = eps;

%% ------------------------------------------------------------ full model
N  = 256;
h  = 1/N;
c  = 1.0;  nu = 1e-3;
dt = 5e-4;                 % forward Euler: c^2*dt <= 2*nu and nu*dt/h^2 <= 1/2 hold
K  = 400;                  % number of time steps (T = 0.2)
xg = (0:N-1)'*h;

I  = eye(N);
Sp = circshift(I, -1, 2);  % (Sp*u)_i = u_{i+1}  (periodic)
Sm = circshift(I,  1, 2);  % (Sm*u)_i = u_{i-1}  (periodic)
A  = -c*(Sp - Sm)/(2*h) + nu*(Sp - 2*I + Sm)/h^2;
M  = I + dt*A;             % forward Euler map

assert(max(abs(eig(M))) <= 1 + 1e-12, 'forward Euler map is unstable');

x0  = double(xg >= 0.2 & xg < 0.4);        % square pulse
x0g = exp(-((xg - 0.7)/0.05).^2);          % smooth IC, used for one control case

% trajectories: pulse (K steps), pulse (2K steps), Gaussian (K steps); columns x_0,...,x_nsteps
ics = {x0, x0, x0g};
ns  = [K, 2*K, K];
trajs = cell(1,3);
for j = 1:3
    Z = zeros(N, ns(j)+1);
    Z(:,1) = ics{j};
    for k = 1:ns(j)
        Z(:,k+1) = M*Z(:,k);
    end
    trajs{j} = Z;
end
Xall = trajs{1};                           % x_0 ... x_K
X    = Xall(:,1:end-1);                    % regression input  x_0..x_{K-1}
Xp   = Xall(:,2:end);                      % regression output x_1..x_K

%% ------------------------------------------------------------ helpers
opinf0 = @(U,Up)     (U'\Up')';                                        % unregularized least squares
opinfR = @(U,Up,lam) ((U*U' + lam*eye(size(U,1))) \ (U*Up'))';         % Tikhonov
relerr = @(Mh,G)     norm(Mh - G,'fro')/norm(G,'fro');

%% ------------------------------------------------------------ bases / data for each case
rng(0);
noise = 1e-8;
Xall_noisy = Xall + noise*randn(size(Xall));
Xn  = Xall_noisy(:,1:end-1);
Xpn = Xall_noisy(:,2:end);

[Phi_same ,~,~] = svd(X,'econ');                              % THEOREM: basis from regression input
[Phi_withK,~,~] = svd(Xall,'econ');                           % basis also sees x_K
[Phi_long ,~,~] = svd(trajs{2}(:,1:end-1),'econ');            % basis from a longer trajectory
[Phi_other,~,~] = svd(trajs{3}(:,1:end-1),'econ');            % basis from a different initial condition
[Phi_rand ,~]   = qr(randn(N,60),0);                          % arbitrary orthonormal subspace
[Phi_noisy,~,~] = svd(Xn,'econ');                             % noisy data, basis from the same noisy data

sigma = svd(X);
lam   = 1e-8*sigma(1)^2;

names = {'same data (theorem)', 'basis also sees x_K', 'basis from 2K-step traj.', ...
         'basis from other IC (Gauss)', 'arbitrary orthonormal subsp.', ...
         'noisy data (1e-8)', 'Tikhonov (lam=1e-8 s1^2)'};
Phis = {Phi_same, Phi_withK, Phi_long, Phi_other, Phi_rand, Phi_noisy, Phi_same};
Xin  = {X,  X,  X,  X,  X,  Xn,  X};
Xout = {Xp, Xp, Xp, Xp, Xp, Xpn, Xp};
lams = [0, 0, 0, 0, 0, 0, lam];

r_list = [2 4 6 8 10 15 20 30];
nr = numel(r_list);  nc = numel(names);
err = zeros(nr, nc);
for i = 1:nr
    r = r_list(i);
    for j = 1:nc
        V = Phis{j}(:,1:r);
        U = V'*Xin{j};  Up = V'*Xout{j};
        if lams(j) == 0
            Mhat = opinf0(U, Up);
        else
            Mhat = opinfR(U, Up, lams(j));
        end
        err(i,j) = relerr(Mhat, V'*M*V);
    end
end

%% ------------------------------------------------------------ extra diagnostics for the theorem case
orth = zeros(nr,1); closure = zeros(nr,1); memory = zeros(nr,1); cond_u = zeros(nr,1); ident = zeros(nr,1);
for i = 1:nr
    r = r_list(i);
    V = Phi_same(:,1:r);
    U = V'*X;  Up = V'*Xp;
    G = V'*M*V;
    % key step of the proof: W*U' = (X - V*U)*U' = 0
    orth(i) = norm((X - V*U)*U','fro') / (norm(X,'fro')*norm(U,'fro'));
    % closure error: projected trajectory vs. trajectory of the intrinsic reduced model
    Xt = zeros(r, K);
    Xt(:,1) = V'*x0;
    for k = 1:K-1
        Xt(:,k+1) = G*Xt(:,k);
    end
    closure(i) = norm(U - Xt,'fro')/norm(U,'fro');
    % size of the memory/orthogonal contribution V'*M*Vperp*W = Up - G*U
    memory(i) = norm(Up - G*U,'fro')/norm(Up,'fro');
    cond_u(i) = sigma(1)/sigma(r);
end

% identity Mhat - G = (Up - G*U)*U'*(U*U')^{-1}, checked for the arbitrary-subspace case
for i = 1:nr
    r = r_list(i);
    V = Phi_rand(:,1:r);
    U = V'*X;  Up = V'*Xp;
    G = V'*M*V;
    Mhat = opinf0(U, Up);
    pol  = (Up - G*U)*pinv(U);              % = (Up - G*U)*U'*(U*U')^{-1}, without squaring cond(U)
    ident(i) = norm((Mhat - G) - pol,'fro') / max(norm(Mhat - G,'fro'), 1e-300);
end

%% ------------------------------------------------------------ report
fprintf('N=%d, K=%d, dt=%g, c=%g, nu=%g, noise=%g, numerical rank of X (tol 1e-14*s1): %d\n', ...
        N, K, dt, c, nu, noise, sum(sigma > 1e-14*sigma(1)));

fprintf('\nRelative error ||Mhat - V''MV||_F / ||V''MV||_F\n');
fprintf('%3s', 'r');
for j = 1:nc, fprintf(' %29s', names{j}); end
fprintf('\n');
for i = 1:nr
    fprintf('%3d', r_list(i));
    for j = 1:nc, fprintf(' %29.3e', err(i,j)); end
    fprintf('\n');
end

fprintf('\nDiagnostics for the theorem case (basis from X)\n');
fprintf('%3s %14s %14s %14s %14s\n', 'r', '||WU''|| rel', 'cond(U)*eps', 'closure err', 'memory term');
for i = 1:nr
    fprintf('%3d %14.3e %14.3e %14.3e %14.3e\n', r_list(i), orth(i), cond_u(i)*EPS, closure(i), memory(i));
end

fprintf('\nIdentity Mhat - G = (Up - G U) U'' (U U'')^-1 (arbitrary subspace), relative discrepancy:\n  ');
for i = 1:nr, fprintf('r=%d: %.1e', r_list(i), ident(i)); if i < nr, fprintf(', '); end, end
fprintf('\n');

%% ------------------------------------------------------------ figure
figure('Position', [100 100 1300 480]);
subplot(1,2,1);
mk = {'o-','s-','^-','v-','d-','x-','p-'};
for j = 1:nc
    semilogy(r_list, max(err(:,j), 1e-17), mk{j}, 'LineWidth', 1.2); hold on;
end
grid on; xlabel('reduced dimension r');
ylabel('||M_{hat} - V^T M V||_F / ||V^T M V||_F');
title('OpInf operator error');
legend(names, 'Location', 'southeast', 'FontSize', 7);

subplot(1,2,2);
semilogy(1:60, sigma(1:60)/sigma(1), 'k.-'); hold on;
semilogy(r_list, closure, 'ro-', 'LineWidth', 1.2);
semilogy(r_list, memory,  'bs-', 'LineWidth', 1.2);
grid on; xlabel('index / reduced dimension r');
title('Non-Markovian effects persist in the theorem case');
legend({'\sigma_i/\sigma_1 of X', 'closure error (theorem case)', 'memory term, relative size'}, 'FontSize', 8);

print(gcf, '-dpng', '-r150', 'advdiff_opinf_experiment_matlab.png');
fprintf('\nsaved advdiff_opinf_experiment_matlab.png\n');