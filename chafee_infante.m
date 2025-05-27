clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

N = 128;
% N = 12;
% N = 64;

%% Chafee-Infante as in https://doi.org/10.1137/19M1292448 combined with 
% setup of Allen-Cahn in https://doi.org/10.1016/j.cma.2022.115836

dt = 1e-5;
% t_end = 4;
t_end = 0.1;
% t_end = 1*dt;
% t_end = 20*dt;
nt = round(t_end/dt);

is = [1 2 3];
A1 = diag(ones(N-1,1),-1) -eye(N);
A1 = (A1+A1')*N^2;
A1 = A1 + eye(N); % additional linear term missing in Allen-Cahn description
% JN3 = power2kron(N,3); % J_N^(3)
% A3 = vecwise_kron(eye(N),3)';
% A3 = A3*JN3;
B = zeros(N,1);
B(1) = 1/dt;
Nu = 1;

% boundary conditions
BC = eye(N);
% x(0,t) = u(t)
A1(1,:) = 0;
A1(1,1) = -1/dt;
% A3(1,:) = 0;
BC(1,:) = 0;

% ddt x(1,t) = 0;
A1(end,:) = 0;
% A3(end,:) = 0;
BC(end,:) = 0;

F1 = @(x1) A1*x1;
F3 = @(x1,x2,x3) -BC*(x1.*x2.*x3);

F1X = @(X) F1(X(:,1));
F3X = @(X) F3(X(:,1),X(:,2),X(:,3));

% IN3 = kron2power(N,3); % I_N^(3)
% f = @(x,u) A1*x + A3*IN3*vecwise_kron(x,3) + B*u; % slow!
f = @(x,u) F1(x) + F3(x,x,x) + B*u;

x0 = zeros(N,1);
u_val = @(t) 10*(sin(pi*t)+1); % U_val
% u_val = @(t) 10*(cos(pi*t)+1); % U_val


%% generate ROM basis construction data
X_b = zeros(N,nt+1);
U_b = zeros(1,nt+1);

t = 0;
x = x0;
u = u_val(t);

X_b(:,1) = x0;
U_b(:,1) = u;

for i=1:nt
    x = x + dt*f(x,u);
    t = t + dt;
    u = u_val(t);

    X_b(:,i+1) = x;
    U_b(:,i+1) = u;
end
%% construct ROM basis via POD
[V,S,~] = svd(X_b,'econ');
n = 14;
Vn = V(:,1:n);

%% construct intrusive operators
tA1 = Vn'*A1*Vn;

n2 = n*(n+1)/2;
tA2 = zeros(n,n2);

Jn3 = power2kron(n,3);
tA3 = precompute_rom_operator(F3X,Vn,3)*Jn3;
tB = Vn'*B;

tO = [tB tA1 tA2 tA3];


%% generate rank-sufficient snapshot data
tX0_pure = rank_suff_basis(n,is);
U0_pure = 1;
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(Nu+1:end,:);
U0 = XU(1:Nu,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

% compute time step estimate (3.10)
dt1 = dt_estimate(X_b,U_b,Vn(:,1),dt,is);

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

ns = 1:n;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A3_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);

n_is__ = n_is(n,is);
offset = cumsum(n_is__);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    % ks = [1:p+n_is_(1), p+n_is__(1)+1:p+n_is__(1)+n_is_(2)];
    ks = 1:Nu+n_is_(1);
    for jj = 2:numel(is)
        ks = [ks Nu+offset(jj-1)+(1:n_is_(jj))];
    end

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
    % hA1 = O(:,A_inds(1,1):A_inds(1,2));
    % hA3 = O(:,A_inds(2,1):A_inds(2,2));
    % hB = O(:,B_inds(1):B_inds(2));
    % 
    % B_errors(j) = norm(tB(1:n_,:) - hB);
    % A1_errors(j) = norm(tA1(1:n_,1:n_is_(1)) - hA1);
    % A3_errors(j) = norm(tA3(1:n_,1:n_is_(2)) - hA3); 

    % tB_ = tB(1:n_,:);
    % tA1_ = tA1(1:n_,1:n_is_(1));
    % tA2_ = tA2(1:n_,1:n_is_(2));
    % tA3_ = tA3(1:n_,1:n_is_(3));
    % 
    % tO_ = [tB_ tA1_ tA2_ tA3_];

    tO_ = tO(1:n_,ks);

    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;
end

figure
hold on
semilogy(ns,A1_errors,'x-', 'LineWidth', 2,'DisplayName',"A_1")
semilogy(ns,A3_errors,'x-', 'LineWidth', 2,'DisplayName',"A_3")
semilogy(ns,B_errors,'x-', 'LineWidth', 2,'DisplayName',"B")
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")

figure
hold on
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O chafee-infante")

ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")

save("data/data_chafee_infante","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

