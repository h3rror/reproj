clear all;
close all;

rng(1); % for reproducibility

warning("eliminate input signal u!")

N = 512;
% N = 12;
% N = 64;

%% ice sheet model in https://icerm.brown.edu/video_archive/3942

% dt = 1e-5;
% dt = 1e-3;
dt = .001;
% t_end = 2;
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
D(1,1) = -1; D(end,end) = 1; % somehow Neumann?
% D(1,end) = -1; D(end,1) = 1; % periodic BC
D = sparse(D);
D = D/(2*1000/N);
% K = vecwise_kron(I)';
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
% fx0 = @(xs) 2.5 + 26*(1-exp(((xs-500)/500/4).^2));
fx0 = @(xs) 1e-2 + 630*(xs/2000+.25).^4.*(xs/2000-.75).^4;

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
% a = -a
% a = 4*a
a = round(a)
% f = @(x,u) a*(f3(x) + f8(x));

% c1 = 8.9271e18;
c1 = 8.9271e-13;
c2 = 2.845713606598e7;
f = @(x,u) c1*f3(x) + c2*f8(x);

% generateFOMdata = true
generateFOMdata = false

if generateFOMdata
    %% generate ROM basis construction data
    X_b = zeros(N,nt);
    U_b = zeros(1,nt);

    t = 0;
    x = x0;
    u = u_val(t);

    X_b(:,1) = x0;
    U_b(:,1) = u;

    for i=2:nt
        % x = x + dt*f(x,u); % forward Euler
        x = fsolve(@(x1) x1 - x - dt*f(x1,u),x); % backward Euler
        t = t + dt;
        u = u_val(t);

        X_b(:,i) = x;
        U_b(:,i) = u;
    end

    save("icesheet_FOMdata","X_b")

else
    %% or load FOM data
    load("icesheet_FOMdata","X_b")
    U_b = zeros(1,nt);
end

%% visualize FOM data

% writerObj = VideoWriter("icesheet",'MPEG-4'); %'Motion JPEG AVI');
% writerObj.FrameRate = 15;
% open(writerObj);
% for i =1:100:nt+1
% plot(xs,X_b(:,i));
% ylim([0 2.6])
% xlabel("spatial dimension")
% ylabel("ice thickness")
% legend("t="+num2str((i-1)*dt))
% frame = getframe(gcf);
% writeVideo(writerObj,frame);
% end
% 
% close(writerObj);

%% construct ROM basis via POD
X_POD = X_b(:,1:2001);

[V,S,~] = svd(X_b,'econ');
% V = eye(N); % botch!
n = 7;
Vn = V(:,1:n);

%% generate rank-sufficient snapshot data

tX0_pure = rank_suff_basis(n,is);
% tX0_pure = rank_suff_basis(n,is)/8;
U0_pure = [];
XU = blkdiag(U0_pure,tX0_pure);
tX0 = XU(p+1:end,:);
U0 = XU(1:p,:);

nf = size(XU,2);
tX1 = zeros(n,nf);

% dt1 = 1
% dt1 = 1e14

Vn1 = Vn(:,1);
X1 = X_POD(:,1:end-1);
X2 = X_POD(:,2:end);
dot_X = (X2-X1)/dt;
U_b1 = U_b(:,1:2000);
dt1 = 1/max(abs(sqrt(sum(Vn1'*dot_X.^2,1))./sqrt(sum(getOpInfMatrix(Vn1'*X1,U_b1,is).^2,1))))


for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

tX0 = int32(full(tX0));
U0 = int32(full(U0));

%% construct intrusive operators
% A3 = D*K*kron(D,K*kron(I,I));
Jn3 = power2kron(n,3);
% IN3 = kron2power(N,3);
% tA3 = c1*Vn'*A3*kron(Vn,kron(Vn,Vn))*Jn3; 

% tA3_2 = precompute_rom_operator(F3X,Vn,3)*Jn3;
tA3 = c1*precompute_rom_operator(F3X,Vn,3)*Jn3;
% 
% norm(tA3_2 - tA3)

Jn8 = power2kron(n,8);
tA8 = c2*precompute_rom_operator(F8X,Vn,8)*Jn8;

%%
ns = 1:n;
% ns = n;
nn = numel(ns);

B_errors = zeros(nn,1);
A3_errors = zeros(nn,1);
A8_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);


condsD = zeros(nn,1);

n_is__ = n_is(n,is);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+p;

    ks = [1:p+n_is_(1), p+n_is__(1)+1:p+n_is__(1)+n_is_(2)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
    hA3 = O(:,A_inds(1,1):A_inds(1,2));
    hA8 = O(:,A_inds(2,1):A_inds(2,2));
    % hB = O(:,B_inds);

    % B_errors(j) = norm(tB(1:n_,:) - hB);
    tA3_ = tA3(1:n_,1:n_is_(1));
    A3_errors(j) = norm(tA3_ - hA3)/norm(tA3_);
    % A3_errors(j) = norm(tA3_ - hA3);
    tA8_ = tA8(1:n_,1:n_is_(2));
    A8_errors(j) = norm(tA8_ - hA8)/norm(tA8_);
    % A8_errors(j) = norm(tA8_ - hA8);

    tO_ = [tA3_ tA8_];
    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;
end

figure
hold on
% semilogy(ns,B_errors,'x-', 'LineWidth', 2,'DisplayName',"B")
% semilogy(ns,A3_errors,'x-', 'LineWidth', 2,'DisplayName',"A_3")
% semilogy(ns,A8_errors,'x-', 'LineWidth', 2,'DisplayName',"A_8")
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O icesheet")

ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")
%%

save("data_icesheet","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

