clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

N = 128;

%% Burgers' model based on https://epubs.siam.org/doi/epdf/10.1137/19M1292448

Omega = [-1 1];
xs = linspace(Omega(1),Omega(2),N);
dx = (Omega(2)-Omega(1))/N;

dt = 1e-4;
t_end = 1;
nt = t_end/dt;

is = [1 2];
I = speye(N);

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

%%
F1X = @(X) F1(X(:,1));
F2X = @(X) F2(X(:,1),X(:,2));

Nu = 0; % input signal dimension

f = @(x,u) F1(x) + F2(x,x);


%% generate ROM basis construction data
X_b = zeros(N,nt+1);
U_b = zeros(Nu,nt+1); 
% X0s = 10*[-sin(pi/2*xs)' sin(3*pi/2*xs)']; % -> make intial condition satisfy BC
x0 = -sin(pi/2*xs)' ; % -> make intial condition satisfy BC

t = 0;
x = x0;
u = U_b(:,1);

X_b(:,1) = x0;
% U_b(:,1) = u;

for i=1:nt
    x = x + dt*f(x,u);
    t = t + dt;
    u = U_b(:,i);

    X_b(:,i+1) = x;
    % U_b(:,i+1) = u;
end

%% construct ROM basis via POD
[V,S,~] = svd(X_b,'econ');
n = 10;
Vn = V(:,1:n);

%% construct intrusive operators
tA1 = precompute_rom_operator(F1X,Vn,1);

Jn2 = power2kron(n,2);
tA2 = precompute_rom_operator(F2X,Vn,2)*Jn2;

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
tX1 = zeros(n,nf);

% compute time step estimate (3.10)
dt1 = dt_estimate(X_b,U_b,Vn(:,1),dt,is);

for i = 1:nf
    tX1(:,i) = Vn'*single_step(Vn*tX0(:,i),U0(:,i),dt1,f);
end

dot_tX = (tX1-tX0)/dt1;

tX0 = int32(full(tX0));
U0 = int32(full(U0));

%%
ns = 1:n;
nn = numel(ns);

B_errors = zeros(nn,1);
A1_errors = zeros(nn,1);
A2_errors = zeros(nn,1);

O_errors = zeros(nn,1);
condsD = zeros(nn,1);

n_is__ = n_is(n,is);

h_energy_error = zeros(nn,1);
t_energy_error = zeros(nn,1);

h_symmetry_error = zeros(nn,1);
t_symmetry_error = zeros(nn,1);

for j = 1:nn
    n_ = ns(j);
    n_is_ = n_is(n_,is);
    nf_ = sum(n_is_)+Nu;

    ks = [1:Nu+n_is_(1), Nu+n_is__(1)+1:Nu+n_is__(1)+n_is_(2)];

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);
    hA1_ = O(:,A_inds(1,1):A_inds(1,2));
    hA2_ = O(:,A_inds(2,1):A_inds(2,2));
    hB_ = O(:,B_inds(1,1):B_inds(1,2));

    % tB_ = tB(1:n_,:);
    tA1_ = tA1(1:n_,1:n_is_(1));
    tA2_ = tA2(1:n_,1:n_is_(2));

    tO_ = [tA1_ tA2_];
    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;


    %% compute energy-preserving constraint violation
    % eq. (19) in https://arxiv.org/pdf/2401.02889
    Jn_3 = power2kron(n_,3);
    In_2 = kron2power(n_,2);
    h_conv = hA2_*In_2;
    h_energy_error(j) = sum(abs(Jn_3'*h_conv(:)));
    % h_energy_error(j) = norm((Jn_3'*h_conv(:)));
    t_conv = tA2_*In_2;
    t_energy_error(j) = sum(abs(Jn_3'*t_conv(:)));
    % t_energy_error(j) = norm((Jn_3'*t_conv(:)));

    %% compute symmetry violation
    h_symmetry_error(j) = norm(hA1_ - hA1_')/norm(hA1_);
    t_symmetry_error(j) = norm(tA1_ - tA1_')/norm(tA1_);

    %% plot eigenvalues of diffusion matrix
    figure(316311)
    hold on
    semilogy(n_*ones(n_,1),-eig(tA1_),'bo', "DisplayName","intrusive", "MarkerSize",10)
    semilogy(n_*ones(n_,1),-eig(hA1_),'rx', "DisplayName","exactOpInf", "MarkerSize",10)
    ylabel("negated eigenvalues","Interpreter","latex", "FontSize",15)
    xlabel("ROM dimension","Interpreter","latex", "FontSize",15)
    set(gca, 'YScale', 'log')
    % legend("show")
    ylim([4 2e4])
    grid on
    legend("intrusive","exactOpInf","Location","northwest","Interpreter","latex", "FontSize",12)

end

savefig("figures/eig_vals.fig")
exportgraphics(gcf,"figures/eig_vals.pdf")

% figure
% hold on
% semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf")
% ylabel("operator error")
% xlabel("ROM dimension")
% set(gca, 'YScale', 'log')
% grid on
% legend("show")


figure
hold on
semilogy(ns,h_energy_error,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf", "MarkerSize",10)
semilogy(ns,t_energy_error,'+:', 'LineWidth', 2,'DisplayName',"intrusive", "MarkerSize",10)
ylabel("energy-preserving constraint violation","Interpreter","latex", "FontSize",15)
xlabel("ROM dimension","Interpreter","latex", "FontSize",15)
set(gca, 'YScale', 'log')
grid on
legend("show","Interpreter","latex", "FontSize",12)
legend("Location","northwest")

savefig("figures/energy_violation.fig")
exportgraphics(gcf,"figures/energy_violation.pdf")

figure
hold on
semilogy(ns,h_symmetry_error,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf", "MarkerSize",10)
semilogy(ns,t_symmetry_error,'+:', 'LineWidth', 2,'DisplayName',"intrusive", "MarkerSize",10)
ylabel("diffusion matrix symmetry violation","Interpreter","latex", "FontSize",15)
xlabel("ROM dimension","Interpreter","latex", "FontSize",15)
set(gca, 'YScale', 'log')
grid on
legend("show","Interpreter","latex", "FontSize",12)
legend("Location","northwest")
ylim([1e-17 1e-15])

savefig("figures/symmetry_violation.fig")
exportgraphics(gcf,"figures/symmetry_violation.pdf")


%% visualize singular values
% figure; semilogy(diag(S),'o-')
% hold on

save("data/data_burgers","O_errors","condsD");


%% FOM solver running for one time step
function x_1 = single_step(x_0,u_0,dt,f)
    x_1 = x_0 + dt*f(x_0,u_0);
end

