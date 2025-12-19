clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

N = 128;
% N = 16;
% N = 6;

    % monolithic = false
    monolithic = true


%% 1D heat equation with temperature-dependent and parameter-dependent (Gaussian) conductivity

Omega = [-1 1];
xis = linspace(Omega(1),Omega(2),N)';
dx = (Omega(2)-Omega(1))/N;

% k0 = @(xis,mu) exp(-(xis-mu).^2);
% syms k0(xis_,mu);
syms xis_ mu_;
k0_(xis_,mu_) = exp(-.5*(xis_-mu_).^2);
k0 = matlabFunction(k0_);
nu = @(x,mu) k0(xis,mu).*x; % thermal conductivity 

x0 = -sin(pi/2*xis) + 1; % -> make intial condition satisfy BC
% x0 = ones(size(xis)) ; % -> make intial condition satisfy BC
mu0 = .3;

figure
hold on
plot(xis,x0, "DisplayName","x_0")
plot(xis,k0(xis,mu0), "DisplayName","k_0(\mu_0)")
plot(xis,nu(x0,mu0), "DisplayName","\nu(x_0,\mu_0)")
legend('show')

%%
dt = 1e-4;
% t_end = 1;
t_end = 100*dt;
% t_end = 10*dt;
nt = t_end/dt;

is = [2];
I = speye(N);

%% 
D = spdiags([-ones(N,1) ones(N,1)], [-1 1],N,N); % first-order central finite difference
D(1,1) = -1; D(end,end) = 1; % homogeneous Neumann BC
D = D/(2*dx);
D1 = D;

% F1 = @(x,theta) 0;
% F1_exact = @(x) D1*(nu(x,k0(mu)).*(D1*x));

F2_exact = @(x1,x2,mu) D1*(nu(x1,mu).*(D1*x2));
F2X_exact = @(X,mu) F2_exact(X(:,1),X(:,2),mu); % enable storing variables in one matrix

% Q: are boundary conditions in D1 correct like this?

% definition for k0 approximations
% F2_k = @(x1,x2,k,mu) D1*((k(xis,mu).*x1).*(D1*x2));
F2_k = @(x1,x2,k) D1*((k.*x1).*(D1*x2));
% F2X_k = @(X,k) F2_k(X(1,:),X(2,:),k,mu0);

% s_max = N;
s_max = 18; % original
% s_max = 6;
% s_max = 1;
% s_max = 1;
% s_max = 30;
% s_max = 21;
mus = linspace(-1,1,s_max);
mus = flip(mus)
% mus = mu0

%% 1) simple setting: fixed mu
% F2 = @(x1,x2,theta) F2_exact(x1,x2,mu0);
% F2X = @(X,theta) F2(X(:,1),X(:,2),theta);  % enable storing variables in one matrix
% f = @(x,u,mu) F2(x,x);
% 
% mus = mu0;
% s= 1;
% qH = 1;
% Theta_H = 1;
% Thetas{1} = Theta_H;

%% 2) general setting: arbitrary qH

k1_ = diff(k0_,xis_);
% k1 = matlabFunction(k1_);
k1 = eval(k1_(xis,mu0));

qH = s_max;
s = qH;
% mus = (0:s-1)+mu0;
% mus = linspace(-1,1,s);
theta_H = @(mu) (mu'-mu0).^(0:s-1);
Theta_H = theta_H(mus(1:s));
Thetas{1} = Theta_H;

% F2 = @(x1,x2,theta) sum(theta'.*[F2_k(x1,x2,k0(xis,mu0)) F2_k(x1,x2,k1)],2);
[F2,kps,k_sum] = F2_taylor_approx(k0_,mu_,qH,F2_k,xis,mu0);
F2X = @(X,theta) F2(X(:,1),X(:,2),theta); % enable storing variables in one matrix

%% some plots
figure
hold on
plot(xis,F2_exact(x0,x0,mu0), "DisplayName","exact")
plot(xis,F2X([x0,x0],theta_H(mu0)),'--',"DisplayName","Taylor approx "+qH)
title("RHS evaluated at x_0 and  \mu_0")
legend("show")

mu1 = -1;
% mu1 = -mu0;
figure
hold on
plot(xis,F2_exact(x0,x0,mu1), "DisplayName","exact")
plot(xis,F2X([x0,x0],theta_H(mu1)),'--',"DisplayName","Taylor approx "+qH)
title("RHS evaluated at x_0 and \mu="+ num2str(mu1))
legend("show")

%%

% f = @(x,u,mu) F2(x,x,theta_H(mu)); % previously used

f = @(x,u,mu) F2_exact(x,x,mu);
% F2_exact = @(x1,x2,mu) D1*(nu(x1,mu).*(D1*x2));


%%

Nu = 0; % input signal dimension

%% generate ROM basis construction data
s_b = 5;

X_b = zeros(N,nt+1,s_b);
U_b = zeros(Nu,nt+1,s_b); 
% X0s = 10*[-sin(pi/2*xs)' sin(3*pi/2*xs)']; % -> make intial condition satisfy BC
% x0 = -sin(pi/2*xis); % -> make intial condition satisfy BC
% mus_b = mus; % so far used
mus_b = linspace(-1,1,s_b);


for k = 1:s_b
    mu = mus_b(:,k);
    X_b(:,:,k) = simulate(x0,dt,nt,@(x) single_step(x,0,dt,f,mu));
end

%% construct ROM basis via POD
[V,S,~] = svd(X_b(:,:),'econ');
n = 30;
% n = 6;
% n = 16;
% n = N;
% n = s_max;
% fac = 12;
% fac = 1;
% n = fac*s_max;

Vn = V(:,1:n);
% Vn = eye(n);

%% singular value decay
figure
semilogy(diag(S)/S(1,1))
title("singular value decay")


%% opinf on ROM basis snapshot data
% tX_b = Vn'*X_b;
tX_b = pagemtimes(Vn',X_b);
tX_b1 = tX_b(:,1:end-1,:);
tX_b2 = tX_b(:,2:end,:);
dot_tX_b = (tX_b2-tX_b1)/dt;

[O,A_inds,B_inds,condD] = p_opinf(dot_tX_b,tX_b1,[],is,Thetas,true);
sota.O = O;
sota.condsD = condD;

%% construct intrusive operators
% tA1s = precompute_rom_operator_param(F1X,Vn,1,qA);
% 
% intr.O = tA1s(:,:);

% Jn2 = power2kron(n,2);
% tA2 = precompute_rom_operator(F2X,Vn,2)*Jn2;

% tA2s = precompute_rom_operator_param(@(X,theta) F2X(X),Vn,2,qH);
tA2s = precompute_rom_operator_param(@(X,theta) F2X(X,theta),Vn,2,qH);
% tA2s = precompute_rom_operator_param(@(X,theta) F2X(X,theta),Vn,2,s_max);

intr.O = tA2s(:,:);
intr.Os = tA2s;

Jn2 = power2kron(n,2);
Omu0 = precompute_rom_operator(@(X) F2X_exact(X,mu0),Vn,2)*Jn2;
Omu1 = precompute_rom_operator(@(X) F2X_exact(X,mu1),Vn,2)*Jn2;

% tA2_2= Vn'*C*kron(Vn,Vn)*Jn2;
% norm(tA2-tA2_2)

%% generate rank-sufficient snapshot data
tX0_pure = rank_suff_basis(n,is);
U0_pure = [];
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(Nu+1:end,:);
U0 = XU(1:Nu,:);

% tX0 = rank_suff_basis(n,is);
% U0 = [];

nf = size(tX0,2);
tX1 = zeros(n,nf,s_max);

% compute time step estimate (3.10)
dt1 = dt_estimate(X_b(:,:,1),U_b(:,:,1),Vn(:,1),dt,is); % internally computes derivatives, so we cannot concatenate trajectories

for k =1:s_max
    mu = mus(:,k);
    for i = 1:nf
        tX1(:,i,k) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f,mu);
    end
end

dot_tX = (tX1-tX0)/dt1;

%%
ns = 1:n;
ss = 1:s_max;
% ns = fac*(1:s_max);
% ns = n;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A2_errors = zeros(nn,1);

deco.O_errors = zeros(nn,1);
deco.condsD = zeros(nn,s);

mono.O_errors = zeros(nn,1);
mono.condsD = zeros(nn,1);

deco.Omu0errors = zeros(nn,1);
deco.Omu1errors = zeros(nn,1);

mono.Omu0errors = zeros(nn,1);
mono.Omu1errors = zeros(nn,1);

%% compute ROM state error
Xmu0_FOM = simulate(x0,dt,nt,@(x) single_step(x,0,dt,f,mu0));
Xmu1_FOM = simulate(x0,dt,nt,@(x) single_step(x,0,dt,f,mu1));

deco.ROMerror_mu0 = zeros(nn,1);
deco.ROMerror_mu1 = zeros(nn,1);

mono.ROMerror_mu0 = zeros(nn,1);
mono.ROMerror_mu1 = zeros(nn,1);

% figure(7)
% hold on
% figure(8)
% hold on
%%

n_is__ = n_is(n,is);

for j = 1:nn
    n_ = ns(j);
    % n_ = 4;
    % n_ = s_max; % botch
    %% NEW: changing s
    % s = n_;
    % s = 1;
    s = qH;
    % s = ns(j); % botch
    % s = ss(j); % botch
    theta_H = @(mu) (mu'-mu0).^(0:s-1);
    Theta_H = theta_H(mus(1:s));
    Thetas{1} = Theta_H;
    %%
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    % ks = [1:Nu+n_is_(1), Nu+n_is__(1)+1:Nu+n_is__(1)+n_is_(2)];
    ks = [1:Nu+n_is_(1)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks,1:s);
    U0_ = U0(:,ks);

    hA2s_ = zeros(n_,n_is_,s);
    tA2s_ = zeros(n_,n_is_,s);

    %% monolithic construction
    if monolithic
        tX = repmat(full(tX0_),1,1,s);
        [O,A_inds,B_inds,condD] = p_opinf(dot_tX_,tX,U0_,is,Thetas,true);
        mono.O = O;
        mono.condsD(j) = condD;
    end

    for k =1:s
        dot_tX_ = dot_tX(1:n_,ks,k);

        [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
        % hA1_ = O(:,A_inds(1,1):A_inds(1,2));
        hA2_ = O(:,A_inds(1,1):A_inds(1,2));
        % hA2_ = O(:,A_inds(2,1):A_inds(2,2));
        hB_ = O(:,B_inds(1,1):B_inds(1,2));

        % tB_ = tB(1:n_,:);
        % tA1_ = tA1s(1:n_,1:n_is_(1),k);
        tA2_ = tA2s(1:n_,1:n_is_(1),k);
        % tA2_ = tA2(1:n_,1:n_is_(2));

        % tO_ = [tA1_ tA2_];
        tO_ = [tA2_];
        % O_errors(j,k) = norm(O-tO_,"fro")/norm(tO_,"fro");

        deco.condsD(j,k) = condD;

        hA2s_(:,:,k) = hA2_;
        tA2s_(:,:,k) = tA2_;
    end

    % % deco.O = hA2s_(:,:)*kron(inv(Theta_H),eye(n_is_))';
    deco.O = hA2s_(:,:)/kron((Theta_H),eye(n_is_))';
    % % O_errors(j) = norm(tA1s_(:,:)-hA1s_,"fro")/norm(tO_,"fro");
    % deco.O_errors(j) = norm(tA2s_(:,:)-deco.O,"fro")/norm(tA2s_(:,:), "fro");
    % if monolithic
    %     mono.O_errors(j) = norm(tA2s_(:,:)-mono.O,"fro")/norm(tA2s_(:,:), "fro");
    % end

    deco.Os = reshape(deco.O,[n_, n_is_(1), s]);

    Omu1_ = Omu1(1:n_,1:n_is_(1));
    deco.Omu1errors(j) = norm(affine_op(deco.Os,theta_H(mu1))-Omu1_,"fro")/norm(Omu1,"fro");
    Omu0_ = Omu0(1:n_,1:n_is_(1));
    deco.Omu0errors(j) = norm(affine_op(deco.Os,theta_H(mu0))-Omu0_,"fro")/norm(Omu0,"fro");
    if monolithic
        mono.Os = reshape(mono.O,[n_, n_is_(1), s]);
        mono.Omu1errors(j) = norm(affine_op(mono.Os,theta_H(mu1))-Omu1_,"fro")/norm(Omu1,"fro");
        mono.Omu0errors(j) = norm(affine_op(mono.Os,theta_H(mu0))-Omu0_,"fro")/norm(Omu0,"fro");
    end

    deco.O_errors(j) = sum(pagenorm(tA2s_-deco.Os,"fro"))/sum(pagenorm(tA2s_, "fro"));
    if monolithic
        mono.O_errors(j) = sum(pagenorm(tA2s_-mono.Os,"fro"))/sum(pagenorm(tA2s_, "fro"));
    end

    %% compute ROM state error
    % Xmu0_FOM = simulate(x0,dt,nt,@(x) single_step(x,0,dt,f,mu0));
    % Xmu1_FOM = simulate(x0,dt,nt,@(x) single_step(x,0,dt,f,mu0));

    In_2 = kron2power(n_,2);
    % ROMmu0_step = @(x) x+ dt*Omu0_*In_2*kron(x,x);
    deco.Omu0_ = affine_op(deco.Os,theta_H(mu0));
    deco.ROMmu0_step = @(x) x+ dt*deco.Omu0_*In_2*kron(x,x);
    % ROMmu1_step = @(x) x+ dt*Omu1_*In_2*kron(x,x);
    deco.Omu1_ = affine_op(deco.Os,theta_H(mu1));
    deco.ROMmu1_step = @(x) x+ dt*deco.Omu1_*In_2*kron(x,x);

    Vn_ = Vn(:,1:n_);
    deco.Xmu0 = simulate(Vn_'*x0,dt,nt,deco.ROMmu0_step);
    deco.Xmu1 = simulate(Vn_'*x0,dt,nt,deco.ROMmu1_step);

    deco.ROMerror_mu0(j) = sum(vecwise_2norm(Xmu0_FOM - Vn_*deco.Xmu0))/sum(vecwise_2norm(Xmu0_FOM));
    deco.ROMerror_mu1(j) = sum(vecwise_2norm(Xmu1_FOM - Vn_*deco.Xmu1))/sum(vecwise_2norm(Xmu1_FOM));
     
    if monolithic
        % ROMmu0_step = @(x) x+ dt*Omu0_*In_2*kron(x,x);
        mono.Omu0_ = affine_op(mono.Os,theta_H(mu0));
        mono.ROMmu0_step = @(x) x+ dt*mono.Omu0_*In_2*kron(x,x);
        % ROMmu1_step = @(x) x+ dt*Omu1_*In_2*kron(x,x);
        mono.Omu1_ = affine_op(mono.Os,theta_H(mu1));
        mono.ROMmu1_step = @(x) x+ dt*mono.Omu1_*In_2*kron(x,x);

        Vn_ = Vn(:,1:n_);
        mono.Xmu0 = simulate(Vn_'*x0,dt,nt,mono.ROMmu0_step);
        mono.Xmu1 = simulate(Vn_'*x0,dt,nt,mono.ROMmu1_step);

        mono.ROMerror_mu0(j) = sum(vecwise_2norm(Xmu0_FOM - Vn_*mono.Xmu0))/sum(vecwise_2norm(Xmu0_FOM));
        mono.ROMerror_mu1(j) = sum(vecwise_2norm(Xmu1_FOM - Vn_*mono.Xmu1))/sum(vecwise_2norm(Xmu1_FOM));
    end

    % ts = 0:dt:t_end;
    % figure(7)
    % semilogy(ts,ROMerror_mu0,'DisplayName',"n=" +num2str(n_))
    % figure(8)
    % semilogy(ts,ROMerror_mu1,'DisplayName',"n=" +num2str(n_))
end

% figure(7)
% ylabel("ROM state error")
% xlabel("time")
% set(gca, 'YScale', 'log')
% grid on
% legend("show")
% title("\mu_0")
% 
% figure(8)
% ylabel("ROM state error")
% xlabel("time")
% set(gca, 'YScale', 'log')
% grid on
% legend("show")
% title("\mu_1")

figure
hold on
semilogy(ns,sum(deco.O_errors,2)/qH,'x-', 'LineWidth', 2,'DisplayName',"decoupled")
if monolithic
    semilogy(ns,sum(mono.O_errors,2)/qH,'+--', 'LineWidth', 2,'DisplayName',"monolithic")
end
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show","Location","southeast")
% title("compares against Taylor approximation!")
exportgraphics(gcf,"figures_deco/taylor_errors.pdf")


figure
hold on
if monolithic
    semilogy(ns,mono.condsD,'x-', 'LineWidth', 2,'DisplayName',"monolithic")
end
semilogy(ns,sum(deco.condsD,2)/s,'x-', 'LineWidth', 2,'DisplayName',"decoupled")
% semilogy(n,sota.condsD, 'x-', 'LineWidth', 2,'DisplayName',"state of the art")
ylabel("condition number")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show","Location","east")
exportgraphics(gcf,"figures_deco/condition_numbers.pdf")


figure
hold on
semilogy(ns,deco.Omu0errors,'x-','DisplayName',"\mu_0 decoupled")
semilogy(ns,deco.Omu1errors,'x-','DisplayName',"\mu_1 decoupled")
if monolithic
    semilogy(ns,mono.Omu0errors,'+--','DisplayName',"\mu_0 monolithic")
    semilogy(ns,mono.Omu1errors,'+--','DisplayName',"\mu_1 monolithic")
end
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show","Location","southeast")
% title("compares against exact Gaussian!")
exportgraphics(gcf,"figures_deco/gaussian_errors.pdf")

figure
hold on
semilogy(ns,deco.ROMerror_mu0,'x-','DisplayName',"\mu_0 decoupled")
semilogy(ns,deco.ROMerror_mu1,'x-','DisplayName',"\mu_1 decoupled")
if monolithic
    semilogy(ns,mono.ROMerror_mu0,'+--','DisplayName',"\mu_0 monolithic")
    semilogy(ns,mono.ROMerror_mu1,'+--','DisplayName',"\mu_1 monolithic")
end
ylabel("average ROM state error")
xlabel("increasing dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")
% title("ROM dim and Taylor dim increasing")
% title("ROM dim increasing, Taylor dim = "+num2str(s))
% title("ROM dim =" + num2str(n_) + ", Taylor dim increasing")
exportgraphics(gcf,"figures_deco/rom_state_errors.pdf")


%% visualize singular values
% figure; semilogy(diag(S),'o-')
% hold on

% save("data/data_burgers","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f,mu)
    x_1 = x_0 + dt*f(x_0,u_0,mu);
end

