clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

% N = 128; % 2^7
N = 2^10;
% N = 2^12;
% N = 12;
% N = 64;

dx = 1/N;

%% 

is = [1];
Nu = 0;
% Nu = 1;

A1_diff = diag(ones(N-1,1),-1) -eye(N);
A1_diff = (A1_diff+A1_diff');

bc_type = "periodic"
% bc_type = "hom_Neumann"
% bc_type = "Dirichlet10"

if bc_type == "hom_Neumann"
    %% homogeneous Neumann BC
    A1_diff(1,1) = -1;
    A1_diff(end,end) = -1;
elseif bc_type == "periodic"
    %% periodic BC
    A1_diff(1,end) = 1;
    A1_diff(end,1) = 1;
    B_diff = 0;
elseif bc_type =="Dirichlet10"
    %% Dirichlet BC: left 1, right 0
    % A1_diff(1,:) = 0; destroys symmetry!
    B_diff = zeros(N,1); B_diff(1)=1;
else
    error("unknown bc_type")
    %%
end
A1_diff = A1_diff/dx^2;
B_diff = B_diff/dx^2;

A1_adv = diag(ones(N-1,1),-1) - diag(ones(N-1,1),1);
if bc_type == "hom_Neumann"
    %% homogeneous Neumann BC
    A1_adv(1,1) = 1;
    A1_adv(end,end) = -1;
elseif bc_type == "periodic"
    %% periodic BC
    A1_adv(1,end) = 1;
    A1_adv(end,1) = -1;
    B_adv = 0;
elseif bc_type =="Dirichlet10"
    %% Dirichlet BC: left 1, right 0
    % A1_adv(1,:) = 0; destroys skew-symmetry!
    B_adv = zeros(N,1); B_adv(1) = 1;
else
    error("unknown bc_type")
    %%
end
% v_adv = 1/dx; % advection velocity: here tuned to be same order of magnitude as diffusion
v_adv = 1; % discretization-independent advection velocity
A1_adv = v_adv*A1_adv/(2*dx);
B_adv = v_adv*B_adv/(2*dx);

% nu = 1;
% nu = .5;
nu = .5*dx;
% nu = 0;
% nu = 0.1;
A1 = nu*A1_diff + A1_adv;
B = nu*B_diff + B_adv;

%% alternative: upwinding
% A1_upw = diag(ones(N-1,1),-1) - eye(N);
% if bc_type == "periodic"
%     %% periodic BC
%     A1_upw(1,end) = 1;
% else
%     error("unknown bc_type")
%     %%
% end
% A1_upw = A1_upw/dx;
%%

% A1 = A1_upw;

% % boundary conditions
% BC = eye(N);
% 
% % x(0,t) = u(t)
% A1(1,:) = 0;
% A1(1,1) = -1/dt;
% BC(1,:) = 0;
% 
% % ddxi x(1,t) = 0
% A1(end,end) = -1/dx^2 + 1;

%% choose time step size

% dt = 1e-5;
% % dt = 1e-2;
% % dt = .5*dx;
% % dt = .25*dx;
% % t_end = 4;
% % t_end = 1;
% % t_end = 2;
% % t_end = 0.1;
% % t_end = 1*dt;
% % t_end = 200*dt;
% % t_end = 1000*dt;
% t_end = 10000*dt;
% nt = round(t_end/dt);

% alternative
% dt_diff = (dx^2 * dy^2) / (2*D*(dx^2 + dy^2));
dt_diff = dx^2/(2*nu)
dt_adv  = 1 / (abs(v_adv)/dx)
dt = 0.5 * min(dt_diff, dt_adv);
% t_end = 1/v_adv;
t_end = 1/v_adv/2;
% t_end = 1/v_adv/20;
T_final = t_end;
nt = ceil(T_final/dt);
dt = T_final/nt;

%%

% F1 = @(x1,u) A1*x1;
if Nu == 0
    F1 = @(x1,u) A1*x1 + B;
elseif Nu == 1
    F1 = @(x1,u) A1*x1 + B*u;
end

F1X = @(X) F1(X(:,1));
% F1X = @(X) F1(X(:,1),u);

f = @(x,u) F1(x,u);

% x0 = zeros(N,1);
xs = (1:N)/N;
% x0 = exp(-(40*(xs-.5)).^2);
x0 = 0*xs;
% x0(1) = 1;
x0(1:(N/2)) = 1;
% x0(1:(N/8)) = 1;
x0 = x0' ;

u_val = @(t) 1+0*t;

%% generate ROM basis construction data
X_b = zeros(N,nt+1);
U_b = zeros(Nu,nt+1);

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


%% state plots
figure; hold on
plot(X_b(:,1))
plot(X_b(:,2))
plot(X_b(:,3))
plot(X_b(:,5))
plot(X_b(:,10))
plot(X_b(:,100))
plot(X_b(:,end))
legend("show")

%% input signal plots
% figure
% ts = linspace(0,t_end,t_end/dt);
% plot(U_b)

%% construct ROM basis via POD
snapshot_stride = 1;
inds = 1:snapshot_stride:size(X_b,2);

% [V,S,~] = svd(X_b(:,inds),'econ');
[V,S,~] = svd(X_b(:,1:end-1),'econ'); % claude suggestion
% n = 14;
n = 28;
% n = 40;
% n = 26;
Vn = V(:,1:n);

%% singular value decay
figure
semilogy(diag(S)/S(1,1))
title("singular value decay")

%% plot POD modes
figure; hold on
for i = 1:n
% for i = n:n
    plot(Vn(:,i))
end
legend("show")

%% prepare standard opinf
tX_b = Vn' *X_b;
tX_b0 = tX_b(:,1:end-1);
tX_b1 = tX_b(:,2:end);
dot_tX_b = (tX_b1-tX_b0)/dt;
% dot_tX_b = (tX_b1-tX_b0-Vn'*B)/dt; %botch!

%% construct intrusive operators
tA1 = Vn'*A1*Vn;

tB = Vn'*B;

% n2 = n*(n+1)/2;
% tA2 = zeros(n,n2);
% 
% Jn3 = power2kron(n,3);
% tA3 = precompute_rom_operator(F3X,Vn,3)*Jn3;
% tB = Vn'*B;
% 
% tO = [tB tA1 tA2 tA3];

if Nu == 0
    tO = tA1;
elseif Nu == 1
    tO = [tB tA1];
end


%% generate rank-sufficient snapshot data
tX0_pure = rank_suff_basis(n,is);

if Nu == 0
    U0_pure = [];
elseif Nu == 1
    U0_pure = [1];
end

XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(Nu+1:end,:);
U0 = XU(1:Nu,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

%% plot initial conditions

figure
hold on
for i = 1:nf
    plot(Vn*tX0(:,i))
end

%%

% compute time step estimate (3.10)
dt1 = dt_estimate(X_b,U_b,Vn(:,1),dt,is);
% dt1 = dt

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

ns = 1:n;
nn = numel(ns);

% B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
% A3_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);

s_O_errors = zeros(nn,1);
s_condsD = zeros(nn,1);

r_O_errors = zeros(nn,1);
r_condsD = zeros(nn,1);

h_ROM_state_error = zeros(nn,1);
t_ROM_state_error = zeros(nn,1);
s_ROM_state_error = zeros(nn,1);
r_ROM_state_error = zeros(nn,1);
rt_ROM_state_error = zeros(nn,1);

selects = {};

best_approx_error_POD = zeros(nn,1);
best_approx_error_mPOD = zeros(nn,1);

n_is__ = n_is(n,is);
offset = cumsum(n_is__);

for j = 1:nn
% for j = nn:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    % ks = [1:p+n_is_(1), p+n_is__(1)+1:p+n_is__(1)+n_is_(2)];
    ks = 1:Nu+n_is_(1);
    for jj = 2:numel(is)
        ks = [ks Nu+offset(jj-1)+(1:n_is_(jj))];
    end

    %% compute exact opinf
    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);

    tO_ = tO(1:n_,ks);

    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;

    %% compute standard opinf
    U_b0 = U_b(:,1:end-1);
    [sO_,~,~,s_condD] = opinf(dot_tX_b(1:n_,:),tX_b0(1:n_,:),U_b0,is,true);

    s_condsD(j) = s_condD;
   
    s_O_errors(j) = norm(sO_-tO_,"fro")/norm(tO_,"fro");
end

figure
hold on
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O exact opinf")
semilogy(ns,s_O_errors,'x-', 'LineWidth', 2,'DisplayName',"O standard opinf")

ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")
box on