clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

% N = 128;
N = 16;
% N = 6;

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
t_end = 1;
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
% Q: are boundary conditions in D1 correct like this?

% definition for k0 approximations
% F2_k = @(x1,x2,k,mu) D1*((k(xis,mu).*x1).*(D1*x2));
F2_k = @(x1,x2,k) D1*((k.*x1).*(D1*x2));
% F2X_k = @(X,k) F2_k(X(1,:),X(2,:),k,mu0);

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

qH = 5;
s = qH;
% mus = (0:s-1)+mu0;
mus = linspace(-1,1,s);
theta_H = @(mu) (mu'-mu0).^(0:s-1);
Theta_H = theta_H(mus);
Thetas{1} = Theta_H;

% F2 = @(x1,x2,theta) sum(theta'.*[F2_k(x1,x2,k0(xis,mu0)) F2_k(x1,x2,k1)],2);
[F2,kps,k_sum] = F2_taylor_approx(k0_,mu_,qH,F2_k,xis,mu0);
F2X = @(X,theta) F2(X(:,1),X(:,2),theta); % enable storing variables in one matrix

% figure
% hold on
% mus2 = xis;
% plot(mus2,k0(1,mus2), "displayname", "exact")
% k1_ = taylor(k0_,mu_,ExpansionPoint = mu0, Order=1);
% plot(mus2,eval(k1_(1,mus2)), "displayname", "Taylor 1")
% k2_ = taylor(k0_,mu_,ExpansionPoint = mu0, Order=2);
% plot(mus2,eval(k2_(1,mus2)), "displayname", "Taylor 2")
% k3_ = taylor(k0_,mu_,ExpansionPoint = mu0, Order=3);
% plot(mus2,eval(k3_(1,mus2)), "displayname", "Taylor 3")
% 
% legend("show")
% title("\xi=1")
% xlabel("\mu")
% ylabel("k_0")

% plot(mus2,k_sum(theta_H(mus2')'))

% f = @(x,u,mu) F2_exact(x,x,mu);
% f = @(x,u,mu) F2(x,x,theta_H(mu));
k_taylor_ = taylor(k0_,mu_,ExpansionPoint = mu0, Order = qH);
% k_taylor = @(xis,mu) matlabFunction((k_taylor_(xis,mu)));
k_taylor = matlabFunction(k_taylor_);
f = @(x,u,mu) F2_k(x,x,k_taylor(xis,mu));

%%

Nu = 0; % input signal dimension

% f = @(x,u,mu) F1(x,theta_A(mu));

% f = @(x,u,mu) F2(x,x);

%% generate ROM basis construction data
X_b = zeros(N,nt+1,s);
U_b = zeros(Nu,nt+1,s); 
% X0s = 10*[-sin(pi/2*xs)' sin(3*pi/2*xs)']; % -> make intial condition satisfy BC
% x0 = -sin(pi/2*xis); % -> make intial condition satisfy BC

for k = 1:s
t = 0;
x = x0;
u = U_b(:,1,s);

X_b(:,1,k) = x0;
% U_b(:,1) = u;

    mu = mus(:,k);
    for i=1:nt
        x = single_step(x,0,dt,f,mu);
        t = t + dt;
        u = U_b(:,i,k);

        X_b(:,i+1,k) = x;
    end
end

%% construct ROM basis via POD
[V,S,~] = svd(X_b(:,:),'econ');
% n = 20;
n = 6;
% n = 16;

Vn = V(:,1:n);
% Vn = eye(n);

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

intr.O = tA2s(:,:);

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
tX1 = zeros(n,nf,s);

% compute time step estimate (3.10)
dt1 = dt_estimate(X_b(:,:,1),U_b(:,:,1),Vn(:,1),dt,is); % internally computes derivatives, so we cannot concatenate trajectories

for k =1:s
    mu = mus(:,k);
    for i = 1:nf
        tX1(:,i,k) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f,mu);
    end
end

dot_tX = (tX1-tX0)/dt1;

%%
ns = 1:n;
% ns = n;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A2_errors = zeros(nn,1);

deco.O_errors = zeros(nn,1);
deco.condsD = zeros(nn,s);

mono.O_errors = zeros(nn,1);
mono.condsD = zeros(nn,1);

n_is__ = n_is(n,is);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    % ks = [1:Nu+n_is_(1), Nu+n_is__(1)+1:Nu+n_is__(1)+n_is_(2)];
    ks = [1:Nu+n_is_(1)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks,:);
    U0_ = U0(:,ks);

    hA2s_ = zeros(n_,n_is_,s);
    tA2s_ = zeros(n_,n_is_,s);

    tX = repmat(full(tX0_),1,1,s);
    [O,A_inds,B_inds,condD] = p_opinf(dot_tX_,tX,U0_,is,Thetas,true);
    mono.O = O;
    mono.condsD(j) = condD;

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

    hA2s_ = hA2s_(:,:)*kron(inv(Theta_H),eye(n_is_))';
    % O_errors(j) = norm(tA1s_(:,:)-hA1s_,"fro")/norm(tO_,"fro");
    deco.O_errors(j) = norm(tA2s_(:,:)-hA2s_,"fro");
    mono.O_errors(j) = norm(tA2s_(:,:)-mono.O,"fro");

end

figure
hold on
semilogy(ns,sum(mono.O_errors,2)/qH,'x-', 'LineWidth', 2,'DisplayName',"monolithic")
semilogy(ns,sum(deco.O_errors,2)/qH,'x-', 'LineWidth', 2,'DisplayName',"decoupled")
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")

figure
hold on
semilogy(ns,mono.condsD,'x-', 'LineWidth', 2,'DisplayName',"monolithic")
semilogy(ns,sum(deco.condsD,2)/s,'x-', 'LineWidth', 2,'DisplayName',"decoupled")
semilogy(n,sota.condsD, 'x-', 'LineWidth', 2,'DisplayName',"state of the art")
ylabel("condition number")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")



%% visualize singular values
% figure; semilogy(diag(S),'o-')
% hold on

% save("data/data_burgers","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f,mu)
    x_1 = x_0 + dt*f(x_0,u_0,mu);
end

