clear all;
close all;

rng(1); % for reproducibility


N = 128;
% N = 12;
% N = 64;
% N = 6;

%% Burgers' model in https://epubs.siam.org/doi/epdf/10.1137/19M1292448

Omega = [-1 1];
xs = linspace(Omega(1),Omega(2),N);
dx = (Omega(2)-Omega(1))/N;

% dt = 1e-5;
dt = 1e-4;
% dt = .05;
t_end = 1;
% t_end = .1;
% t_end = 1*dt;
% t_end = 20*dt;
nt = t_end/dt;

is = [1 2];
I = speye(N);
% % % % D = -diag(ones(N-1,1),-1) + I; % not sparse!
% % % D = spdiags([-ones(N,1) ones(N,1)], [-1 1],N,N);
% % % D(1,end) = -1; % periodic BC?
% % % D(end,1) = 1; % periodic BC?
% % % % D([1 N],:) = 0; % dx/dt = 0 BC / Dirichlet BC
% % % % D(1,end) = -1; D(end,1) = 1;
% % % 
% % % % D = sparse(D);
% % % D = D/(2*dx);

%% skew-symmetric convection operator for Burger's equation as described in
% https://www.sciencedirect.com/science/article/pii/S0021999124002523
index = @(x) mod(x-1,N)+1;
kron_ind = @(i,j) i+(j-1)*N;

C = zeros(N,N^2);
for i = 1:N
    C(i,kron_ind(index(i-1),index(i-1))) = 1; % u_{i-1}^2
    C(i,kron_ind(index(i+1),index(i+1))) = -1; % u_{i+1}^2
    C(i,kron_ind(index(i-1),i)) = 1; % u_i u_{i-1}
    C(i,kron_ind(index(i+1),i)) = -1; % u_i u_{i+1}
end

C = C/(3*dx);

F2 = @(x1,x2) C*kron(x1,x2);
%% negative semi-definite diffusion operator

D = -2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
D(N,1) = 1;
D(1,N) = 1;
D = D/dx^2;

F1 = @(x) D*x;

% D2 = spdiags([ones(N,1) -2*ones(N,1) ones(N,1)], [-1 0 1],N,N);
% D2([1 N],:) = 0; % dx/dt = 0 BC / Dirichlet BC
% D2 = D2/dx^2;
% mu = .5;
% D2 = mu*D2;
% D2(1,1) = -1/dt; % dx/dt = 0 BC / Dirichlet BC
% D2(end,end) = -1/dt; % dx/dt = 0 BC / Dirichlet BC


% BC = eye(N);
% BC([1 N],:) = 0; % dx/dt = 0 BC / Dirichlet BC
% 
B = zeros(N,1);
% B(1) = 1;
% B(end) = -1;
% B = B/dt;
% 
% F1 = @(x1) D2*x1;
% F2 = @(x1,x2) -x1.*(D*x2);
% 
% K = vecwise_kron(I)';
% A2 = -K*kron(I,D);

% A3 = D*K*kron(D,K*kron(I,I));
% A8 = D*K*kron(D,K*kron(D,K*kron(D,K*kron(I,K*kron(I,K*kron(I,K*kron(I,I)))))));
% F3 = @(x1,x2,x3) D*((D*x1).*x2.*x3);
% F8 = @(x1,x2,x3,x4,x5,x6,x7,x8) ...
%         D*((D*x1).*(D*x2).*(D*x3).*x4.*x5.*x6.*x7.*x8);

F1X = @(X) F1(X(:,1));
F2X = @(X) F2(X(:,1),X(:,2));
% F3X = @(X) F3(X(:,1),X(:,2),X(:,3));
% F8X = @(X) F8(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),X(:,8));

% B = 0;
% p = 1;
p = 0;

%%
% x0 = zeros(N,1);

%%
% a = 12.5; % choose such that CFL condition a <= dx^2/dt/u_max^5 is met!
% f = @(x,u) F1(x) + F2(x,x) + B*u;
f = @(x,u) F1(x) + F2(x,x);


%% generate ROM basis construction data
Mt = 2; % make test and training data differ
% Mt = 1; % make test and training data equal

X_b = zeros(N,nt,Mt);
% U_b = zeros(1,nt);
% U_b = 10*rand(1,nt,Mt);
% U_test = 10*rand(1,nt);
% X0s = [x0, rand(N,Mt-1)];
% X0s = rand(N,Mt);
% X0s = [zeros(N,1) sin(2*pi*xs)'];

U_b = 10*ones(1,nt,Mt); % --> insight: input should not be time-dependent!
X0s = 10*[-sin(pi/2*xs)' sin(3*pi/2*xs)']; % -> make intial condition satisfy BC

for k=1:Mt
    x0 = X0s(:,k);

    t = 0;
    x = x0;
    u = U_b(:,1,k);

    X_b(:,1,k) = x0;
    % U_b(:,1) = u;

    for i=2:nt
        x = x + dt*f(x,u);
        t = t + dt;
        u = U_b(:,i,k);

        X_b(:,i,k) = x;
        % U_b(:,i) = u;
    end

end
%% construct ROM basis via POD
[V,S,~] = svd(X_b(:,:),'econ');
% V = eye(N); % botch!
n = 10;
Vn = V(:,1:n);

%% construct intrusive operators
% A3 = D*K*kron(D,K*kron(I,I));
% Jn3 = power2kron(n,3);
% IN3 = kron2power(N,3);
% tA3 = a*Vn'*A3*kron(Vn,kron(Vn,Vn))*Jn3;

% tA3_2 = precompute_rom_operator(F3X,Vn,3)*Jn3;
% 
% norm(tA3_2 - tA3)

tA1 = precompute_rom_operator(F1X,Vn,1);
Jn2 = power2kron(n,2);
tA2 = precompute_rom_operator(F2X,Vn,2)*Jn2;

% tA2_2= Vn'*K*kron(-eye(N),D)*kron(Vn,Vn)*Jn2;
tA2_2= Vn'*C*kron(Vn,Vn)*Jn2;
norm(tA2-tA2_2)

% Jn8 = power2kron(n,8);
% tA8 = a*precompute_rom_operator(F8X,Vn,8)*Jn8;

tB = Vn'*B;


%% generate rank-sufficient snapshot data

tX0_pure = rank_suff_basis(n,is);
% U0_pure = [1];
U0_pure = [];
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(p+1:end,:);
U0 = XU(1:p,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

% dt1 = 1
Vn1 = Vn(:,1);
X1 = X_b(:,1:end-1);
X2 = X_b(:,2:end);
dot_X = (X2-X1)/dt;
U_b1 = U_b(:,1:end-1);
dt1 = 1/max(abs(sqrt(sum(Vn1'*dot_X.^2,1))./sqrt(sum(getOpInfMatrix(Vn1'*X1,U_b1,is).^2,1))))

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

tX0 = int32(full(tX0));
U0 = int32(full(U0));

%%
ns = 1:n;
% ns = 1:6;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A2_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);

sB_errors = zeros(nn,1);
sA1_errors = zeros(nn,1);
sA2_errors = zeros(nn,1);

n_is__ = n_is(n,is);

tX_train = zeros(nn,nt,nn);
hX_train = zeros(nn,nt,nn);
sX_train = zeros(nn,nt,nn);

train_error_t = zeros(nn,1);
train_error_h = zeros(nn,1);
train_error_s = zeros(nn,1);

tX_test = zeros(nn,nt,nn);
hX_test = zeros(nn,nt,nn);
sX_test = zeros(nn,nt,nn);

test_error_t = zeros(nn,1);
test_error_h = zeros(nn,1);
test_error_s = zeros(nn,1);

h_energy_error = zeros(nn,1);
t_energy_error = zeros(nn,1);

h_symmetry_error = zeros(nn,1);
t_symmetry_error = zeros(nn,1);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+p;

    ks = [1:p+n_is_(1), p+n_is__(1)+1:p+n_is__(1)+n_is_(2)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    % [O,A_inds,B_inds] = opinf(dot_tX_,tX0_,U0_,is);
    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
    hA1_ = O(:,A_inds(1,1):A_inds(1,2));
    hA2_ = O(:,A_inds(2,1):A_inds(2,2));
    hB_ = O(:,B_inds(1,1):B_inds(1,2));

    tB_ = tB(1:n_,:);
    tA1_ = tA1(1:n_,1:n_is_(1));
    tA2_ = tA2(1:n_,1:n_is_(2));

    B_errors(j) = norm(tB(1:n_,:) - hB_);
    A1_errors(j) = norm(tA1(1:n_,1:n_is_(1)) - hA1_);
    A2_errors(j) = norm(tA2(1:n_,1:n_is_(2)) - hA2_);

    % tO_ = [tB_ tA1_ tA2_];
    tO_ = [tA1_ tA2_];
    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;


    %% compute energy-preserving constraint violation
    % eq. (19) in https://arc.aiaa.org/doi/epdf/10.2514/6.2024-1012
    Jn_3 = power2kron(n_,3);
    In_2 = kron2power(n_,2);
    h_conv = hA2_*In_2;
    h_energy_error(j) = norm(Jn_3'*h_conv(:));
    % h_energy_error(j) = norm(Jn_3'*h_conv(:))/norm(hA2_); % relative
    % h_energy_error(j) = max(abs(Jn_3'*h_conv(:))); 
    t_conv = tA2_*In_2;
    t_energy_error(j) = norm(Jn_3'*t_conv(:));
    % t_energy_error(j) = norm(Jn_3'*t_conv(:))/norm(tA2_); % relative
    % t_energy_error(j) = max(abs(Jn_3'*t_conv(:))); 

    %% compute symmetry violation
    % h_symmetry_error(j) = norm(hA1_ - hA1_');
    % t_symmetry_error(j) = norm(tA1 - tA1');

    h_symmetry_error(j) = norm(hA1_ - hA1_')/norm(hA1_);
    t_symmetry_error(j) = norm(tA1 - tA1')/norm(tA1);

    % h_symmetry_error(j) = max(abs(hA1_ - hA1_'),[],"all");
    % t_symmetry_error(j) = max(abs(tA1_ - tA1_'),[],"all");

    % plot eigenvalues of diffusion matrix
    figure(316311)
    hold on
    semilogy(n_*ones(n_,1),-eig(tA1_),'bo', "DisplayName","intrusive")
    semilogy(n_*ones(n_,1),-eig(hA1_),'rx', "DisplayName","exactOpInf")
    ylabel("- eigenvalues")
    xlabel("ROM dimension")
    set(gca, 'YScale', 'log')
    % legend("show")
    ylim([4 2e4])
    grid on
    legend("intrusive","exactOpInf","Location","northwest")


    %% compare test and training errors with standard opinf

    %% perform standard opinf only on the first slice of snapshot data
    Vn_ = V(:,1:n_);
    tX_b = Vn_'*X_b(:,:,1);
    dot_tX_b = (tX_b(:,2:end)-tX_b(:,1:end-1))/dt;

    U_opinf = U_b(:,:,1);
    
    [sO,A_inds,B_inds] = opinf(dot_tX_b,tX_b(:,1:end-1),U_opinf(:,1:end-1),is);
    sA1 = sO(:,A_inds(1,1):A_inds(1,2));
    sA2 = sO(:,A_inds(2,1):A_inds(2,2));
    sB = sO(:,B_inds(1,1):B_inds(1,2));

    sB_errors(j) = norm(tB(1:n_,:) - sB);
    sA1_errors(j) = norm(tA1(1:n_,1:n_is_(1)) - sA1);
    sA2_errors(j) = norm(tA2(1:n_,1:n_is_(2)) - sA2);

    %% training error
    % t = 0;
    % 
    % x0_train = X0s(:,1);
    % tx = Vn_'*x0_train;
    % hx = Vn_'*x0_train;
    % sx = Vn_'*x0_train;
    % 
    % U_train = U_b(:,:,1);
    % X_train = X_b(:,:,1);
    % 
    % u = U_train(:,1);
    % 
    % tX_train(1:n_,1,j) = tx;
    % hX_train(1:n_,1,j) = hx;
    % sX_train(1:n_,1,j) = sx;
    % 
    % % In_2 = kron2power(n_,2);
    % 
    % for i=2:nt
    %     tx = single_step(tx,u,dt,@(x,u) tA1_*x + tA2_*In_2*kron(tx,tx) + tB_*u);
    %     hx = single_step(hx,u,dt,@(x,u) hA1_*x + hA2_*In_2*kron(hx,hx) + hB_*u);
    %     sx = single_step(sx,u,dt,@(x,u) sA1*x + sA2*In_2*kron(sx,sx) + sB*u);
    %     t = t + dt;
    %     u = U_train(:,i);
    % 
    %     tX_train(1:n_,i,j) = tx;
    %     hX_train(1:n_,i,j) = hx;
    %     sX_train(1:n_,i,j) = sx;
    % end
    % 
    % train_error_t(j) = norm(Vn_*tX_train(1:n_,:,j) - X_train,"fro")/norm(X_b,"fro");
    % train_error_h(j) = norm(Vn_*hX_train(1:n_,:,j) - X_train,"fro")/norm(X_b,"fro");
    % train_error_s(j) = norm(Vn_*sX_train(1:n_,:,j) - X_train,"fro")/norm(X_b,"fro");
    % 
    % %% test error
    % t = 0;
    % 
    % % x0_test = ones(N,1);
    % x0_test = X0s(:,Mt);
    % tx = Vn_'*x0_test;
    % hx = Vn_'*x0_test;
    % sx = Vn_'*x0_test;
    % 
    % U_test = U_b(:,:,Mt);
    % u = U_test(:,1);
    % 
    % X_test = X_b(:,:,Mt);
    % 
    % tX_test(1:n_,1,j) = tx;
    % hX_test(1:n_,1,j) = hx;
    % sX_test(1:n_,1,j) = sx;
    % 
    % In_2 = kron2power(n_,2);
    % 
    % for i=2:nt
    %     tx = single_step(tx,u,dt,@(x,u) tA1_*x + tA2_*In_2*kron(tx,tx) + tB_*u);
    %     hx = single_step(hx,u,dt,@(x,u) hA1_*x + hA2_*In_2*kron(hx,hx) + hB_*u);
    %     sx = single_step(sx,u,dt,@(x,u) sA1*x + sA2*In_2*kron(sx,sx) + sB*u);
    %     t = t + dt;
    %     u = U_test(:,i);
    % 
    %     tX_test(1:n_,i,j) = tx;
    %     hX_test(1:n_,i,j) = hx;
    %     sX_test(1:n_,i,j) = sx;
    % end
    % 
    % test_error_t(j) = norm(Vn_*tX_test(1:n_,:,j) - X_test,"fro")/norm(X_b,"fro");
    % test_error_h(j) = norm(Vn_*hX_test(1:n_,:,j) - X_test,"fro")/norm(X_b,"fro");
    % test_error_s(j) = norm(Vn_*sX_test(1:n_,:,j) - X_test,"fro")/norm(X_b,"fro");

end

figure
hold on
semilogy(ns,A1_errors,'x-', 'LineWidth', 2,'DisplayName',"A_1 rank-suff")
semilogy(ns,A2_errors,'x-', 'LineWidth', 2,'DisplayName',"A_2 rank-suff")
semilogy(ns,B_errors,'x-', 'LineWidth', 2,'DisplayName',"B rank-suff")

semilogy(ns,sA1_errors,'x-', 'LineWidth', 2,'DisplayName',"A_1 standard opinf")
semilogy(ns,sA2_errors,'x-', 'LineWidth', 2,'DisplayName',"A_2 standard opinf")
semilogy(ns,sB_errors,'x-', 'LineWidth', 2,'DisplayName',"B standard opinf")
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")

figure
hold on
semilogy(ns,train_error_t,'x-', 'LineWidth', 2,'DisplayName',"intrusive")
semilogy(ns,train_error_h,'o-', 'LineWidth', 2,'DisplayName',"rank-suff")
semilogy(ns,train_error_s,'s-', 'LineWidth', 2,'DisplayName',"standard")
ylabel("average relative state error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
title("training error")
grid on
legend("show")

figure
hold on
semilogy(ns,test_error_t,'x-', 'LineWidth', 2,'DisplayName',"intrusive")
semilogy(ns,test_error_h,'o-', 'LineWidth', 2,'DisplayName',"rank-suff")
semilogy(ns,test_error_s,'s-', 'LineWidth', 2,'DisplayName',"standard")
ylabel("average relative state error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
title("test error")
grid on
legend("show")

figure
hold on
semilogy(ns,h_energy_error,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf")
semilogy(ns,t_energy_error,'+:', 'LineWidth', 2,'DisplayName',"intrusive")
ylabel("energy-preserving constraint violation")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")
legend("Location","northwest")


figure
hold on
semilogy(ns,h_symmetry_error,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf")
semilogy(ns,t_symmetry_error,'+:', 'LineWidth', 2,'DisplayName',"intrusive")
ylabel("diffusion matrix symmetry violation")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")
legend("Location","northwest")
ylim([5e-17 1e-15])


figure
hold on
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O burgers")


%% visualize singular values
figure; semilogy(diag(S),'o-')
hold on

save("data_burgers_periodic","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

