clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

N = 512;
% N = 12;
% N = 64;

%% ice sheet model in https://icerm.brown.edu/video_archive/3942
Omega = [0 1000];
dx = (Omega(2)-Omega(1))/N;

dt = .001;
t_end = 2;
nt = t_end/dt;

is = [3 8];
I = speye(N);
D = spdiags([-ones(N,1) ones(N,1)], [-1 1],N,N); % first-order central finite difference
D(1,1) = -1; D(end,end) = 1; % homogeneous Neumann BC
D = D/(2*dx);

F3 = @(x1,x2,x3) D*((D*x1).*x2.*x3);
F8 = @(x1,x2,x3,x4,x5,x6,x7,x8) ...
        D*((D*x1).*(D*x2).*(D*x3).*x4.*x5.*x6.*x7.*x8);
F3X = @(X) F3(X(:,1),X(:,2),X(:,3));
F8X = @(X) F8(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),X(:,8));

B = 0;
Nu = 0; % input signal dimension


xs = linspace(Omega(1),Omega(2),N)';
fx0 = @(xs) 1e-2 + 630*(xs/2000+.25).^4.*(xs/2000-.75).^4;

x0 = fx0(xs);
u_val = @(t) 0; % U_val

%%
f3 = @(x) D*((D*x).*x.*x);
f8 = @(x) D*((D*x).^3.*x.^5);


% coefficients
rho = 910;
g = 9.81;
beta = 1e16;
gamma = 1e-4;

c1 = rho*g/beta;          % 8.9271e-13
c2 = 2*gamma*rho^3*g^3/5; % 2.845713606598e7

f = @(x,u) c1*f3(x) + c2*f8(x);

% generatePODdata = true 
generatePODdata = false
if generatePODdata

    %% generate ROM basis construction data
    X_b = zeros(N,nt+1);
    U_b = zeros(Nu,nt+1);

    t = 0;
    x = x0;
    u = u_val(t);

    X_b(:,1) = x0;
    U_b(:,1) = u;

    for i=1:nt
        x = fsolve(@(x1) x1 - x - dt*f(x1,u),x); % backward Euler
        t = t + dt;
        u = u_val(t);

        X_b(:,i+1) = x;
        U_b(:,i+1) = u;
    end

    save("data/icesheet_FOMdata","X_b")

else
    %% or load FOM data
    load("data/icesheet_FOMdata","X_b")
    U_b = zeros(1,nt+1);
end

%% visualize FOM data

writerObj = VideoWriter("figures/icesheet",'MPEG-4'); %'Motion JPEG AVI');
writerObj.FrameRate = 15;
open(writerObj);
for i =1:100:(nt+1)
plot(xs,X_b(:,i));
ylim([0 2.6])
xlabel("spatial dimension")
ylabel("ice thickness")
legend("t="+num2str((i-1)*dt))
frame = getframe(gcf);
writeVideo(writerObj,frame);
end

close(writerObj);

%% construct ROM basis via POD
X_POD = X_b(:,1:2001);

[V,S,~] = svd(X_b,'econ');
n = 7;
Vn = V(:,1:n);

%% generate rank-sufficient snapshot data

tX0_pure = rank_suff_basis(n,is);
U0_pure = [];
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(Nu+1:end,:);
U0 = XU(1:Nu,:);

nf = size(XU,2);
tX1 = zeros(n,nf);


% compute time step estimate (3.10)
dt1 = dt_estimate(X_b,U_b,Vn(:,1),dt,is);

%%

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

tX0 = int32(full(tX0));
U0 = int32(full(U0));

%% construct intrusive operators
Jn3 = power2kron(n,3);
tA3 = c1*precompute_rom_operator(F3X,Vn,3)*Jn3;

Jn8 = power2kron(n,8);
tA8 = c2*precompute_rom_operator(F8X,Vn,8)*Jn8;

%%
ns = 1:n;
nn = numel(ns);

A3_errors = zeros(nn,1);
A8_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);

n_is__ = n_is(n,is);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    ks = [1:Nu+n_is_(1), Nu+n_is__(1)+1:Nu+n_is__(1)+n_is_(2)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
    hA3 = O(:,A_inds(1,1):A_inds(1,2));
    hA8 = O(:,A_inds(2,1):A_inds(2,2));

    tA3_ = tA3(1:n_,1:n_is_(1));
    A3_errors(j) = norm(tA3_ - hA3)/norm(tA3_);
    tA8_ = tA8(1:n_,1:n_is_(2));
    A8_errors(j) = norm(tA8_ - hA8)/norm(tA8_);

    tO_ = [tA3_ tA8_];
    if n_ == 1
        disp("for ROM dimension n=1, ||\tilde O\|_2="+norm(tO_))
    end
    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;
end

figure
hold on
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O icesheet")

ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")
%%

save("data/data_icesheet","O_errors","condsD");

%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

