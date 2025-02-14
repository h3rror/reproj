clear all;
close all;

rng(1); % for reproducibility

N = 128;
% N = 12;
% N = 64;

%% ice sheet model in https://icerm.brown.edu/video_archive/3942

% dt = 1e-5;
% dt = 1e-3;
dt = .05;
t_end = 10;
% t_end = .1;
% t_end = 1*dt;
% t_end = 20*dt;
nt = t_end/dt;

is = [3 8];
I = speye(N);
% D = -diag(ones(N-1,1),-1) + I; % not sparse!
D = spdiags([-ones(N,1) ones(N,1)], [-1 1],N,N);
% D(1,end) = -1; % periodic BC?
% D([1 N],:) = 0; % dx/dt = 0 BC
D(1,end) = -1; D(end,1) = 1;
D = sparse(D);
D = D/(2*1000/N);
K = vecwise_kron(I)';
% A3 = D*K*kron(D,K*kron(I,I));
% A8 = D*K*kron(D,K*kron(D,K*kron(D,K*kron(I,K*kron(I,K*kron(I,K*kron(I,I)))))));
F3 = @(x1,x2,x3) D*((D*x1).*x2.*x3);
F8 = @(x1,x2,x3,x4,x5,x6,x7,x8) ...
        D*((D*x1).*(D*x2).*(D*x3).*x4.*x5.*x6.*x7.*x8);
F3X = @(X) F3(X(:,1),X(:,2),X(:,3));
F8X = @(X) F8(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),X(:,8));

B = 0;
p = 0;

%% boundary conditions
% % x(0,t) = u(t)
% A1(1,:) = 0;
% A1(1,1) = -1/dt;
% A3(1,:) = 0;
% 
% % ddt x(1,t) = 0;
% A1(end,:) = 0;
% A3(end,:) = 0;

%%
xs = linspace(0,1000,N)';
fx0 = @(xs) 2.5 + 26*(1-exp(((xs-500)/500/4).^2));

x0 = fx0(xs);
u_val = @(t) 0; % U_val

%%
f3 = @(x) D*((D*x).*x.*x);
f8 = @(x) D*((D*x).^3.*x.^5);

% a = 12.5; % choose such that CFL condition a <= dx^2/dt/u_max^5 is met!
dx = 1000/N;
x_max = max(abs(x0))
dx_max = max(abs(D*x0))
a = dx^2/dt/(x_max^2 + x_max^5*dx_max^2)
f = @(x,u) a*(f3(x) + f8(x));


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
[V,S,~] = svd(X_b,'econ');
V = eye(N); % botch!
n = 6;
Vn = V(:,1:n);

%% construct intrusive operators
A3 = D*K*kron(D,K*kron(I,I));
Jn3 = power2kron(n,3);
% IN3 = kron2power(N,3);
tA3 = Vn'*A3*kron(Vn,kron(Vn,Vn))*Jn3;

% tA3_2 = precompute_rom_operator(F3X,Vn,3)*Jn3;
% 
% norm(tA3_2 - tA3)

Jn8 = power2kron(n,8);
tA8 = precompute_rom_operator(F8X,Vn,8)*Jn8;


%% generate rank-sufficient snapshot data

tX0_pure = rank_suff_basis(n,is);
U0_pure = [];
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(p+1:end,:);
U0 = XU(1:p,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt,f);
end

dot_tX = (tX1-tX0)/dt;

tX0 = int32(full(tX0));
U0 = int32(full(U0));

[O,A_inds,B_inds] = opinf(dot_tX,tX0,U0,is);
hA3 = O(:,A_inds(1,1):A_inds(1,2));
hA8 = O(:,A_inds(2,1):A_inds(2,2));
% hB = O(:,B_inds);

norm(tB - hB)
norm(tA3 - hA3)
norm(tA8 - hA8)

%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

