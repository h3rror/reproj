clear all;
close all;

rng(1); % for reproducibility

% N = 128;
N = 12;
% N = 64;

%% Allen-Cahn equation as in https://doi.org/10.1016/j.cma.2022.115836

dt = 1e-5;
% t_end = .1;
% t_end = 1*dt;
t_end = 20*dt;
nt = t_end/dt;

is = [1 3];
A1 = diag(ones(N-1,1),-1) -eye(N);
A1 = (A1+A1')*N^2;
JN3 = power2kron(N,3); % J_N^(3)
A3 = vecwise_kron(eye(N),3)';
A3 = A3*JN3;
B = zeros(N,1);
B(1) = 1/dt;
p = 1;

% boundary conditions
% x(0,t) = u(t)
A1(1,:) = 0;
A1(1,1) = -1/dt;
A3(1,:) = 0;

% ddt x(1,t) = 0;
A1(end,:) = 0;
A3(end,:) = 0;

IN3 = kron2power(N,3); % I_N^(3)
f = @(x,u) A1*x + A3*IN3*vecwise_kron(x,3) + B*u;

x0 = zeros(N,1);
u_val = @(t) 10*(sin(pi*t)+1); % U_val


%% generate ROM basis construction data
X_b = zeros(N,nt);
U_b = zeros(1,nt);

t = 0;
x = x0;
u = u_val(t);

X_b(:,1) = x0;
U_b(:,1) = u;

for i=2:nt
    x = x + dt*f(x,u);
    t = t + dt;
    u = u_val(t);

    X_b(:,i) = x;
    U_b(:,i) = u;
end
%% construct ROM basis via POD
[V,~,~] = svd(X_b,'econ');
n = 10;
Vn = V(:,1:n);

%% construct intrusive operators
tA1 = Vn'*A1*Vn;
Jn3 = power2kron(n,3);
% In3 = kron2power(n,3);
tA3 = Vn'*A3*IN3*kron(Vn,kron(Vn,Vn))*Jn3;
tB = Vn'*B;

%% generate rank-sufficient snapshot data

tX0_pure = rank_suff_basis(n,is);
U0_pure = 1;
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(p+1:end,:);
U0 = XU(1:p,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt,f);
end

dot_tX = (tX1-tX0)/dt;

ns = 1:10;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A3_errors = zeros(nn,1);

n_is__ = n_is(n,is);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+p;

    ks = [1:p+n_is_(1), p+n_is__(1)+1:p+n_is__(1)+n_is_(2)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds] = opinf(dot_tX_,tX0_,U0_,is);
    hA1 = O(:,A_inds(1,1):A_inds(1,2));
    hA3 = O(:,A_inds(2,1):A_inds(2,2));
    hB = O(:,B_inds);

    B_errors(j) = norm(tB(1:n_,:) - hB);
    A1_errors(j) = norm(tA1(1:n_,1:n_is_(1)) - hA1);
    A3_errors(j) = norm(tA3(1:n_,1:n_is_(2)) - hA3);
end

figure
hold on
semilogy(ns,B_errors,'x-', 'LineWidth', 2,'DisplayName',"B")
semilogy(ns,A1_errors,'x-', 'LineWidth', 2,'DisplayName',"A_1")
semilogy(ns,A3_errors,'x-', 'LineWidth', 2,'DisplayName',"A_3")
ylabel("operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")

%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

