clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

%% 

is = [1];
Nu = 0;

%% insert script
assemble_A_advection_diffusion_2d;
A1 = A;
x0 = C0_vec;
N = Nx*Ny;
nt = Nt;


F1 = @(x1) A1*x1;

F1X = @(X) F1(X(:,1));

f = @(x,u) F1(x);

% x0 = zeros(N,1);
xs = (1:N)/N;
x0 = exp(-(40*(xs-.5)).^2);
x0 = x0' ;

u_val = @(t) [];

%% generate ROM basis construction data
[X_b,U_b] = gen_FOM_data(x0,u_val,f,nt,N,dt);


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

[V,S,~] = svd(X_b(:,inds),'econ');
% n = 14;
% n = 24;
n = 40;
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

%% construct intrusive operators
tA1 = Vn'*A1*Vn;

% n2 = n*(n+1)/2;
% tA2 = zeros(n,n2);
% 
% Jn3 = power2kron(n,3);
% tA3 = precompute_rom_operator(F3X,Vn,3)*Jn3;
% tB = Vn'*B;
% 
% tO = [tB tA1 tA2 tA3];
tO = tA1;


%% generate rank-sufficient snapshot data
tX0_pure = rank_suff_basis(n,is);
U0_pure = [];
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

    tX0_ = tX0(1:n_,ks);
    dot_tX_ = dot_tX(1:n_,ks);
    U0_ = U0(:,ks);

    [O,A_inds,B_inds,condD] = opinf(dot_tX_,tX0_,U0_,is,true);

    tO_ = tO(1:n_,ks);

    O_errors(j) = norm(O-tO_,"fro")/norm(tO_,"fro");

    condsD(j) = condD;

    %% compute standard opinf
    [sO_,~,~,s_condD] = opinf(dot_tX_b(1:n_,:),tX_b0(1:n_,:),U_b,is,true);

    s_condsD(j) = s_condD;
   
    s_O_errors(j) = norm(sO_-tO_,"fro")/norm(tO_,"fro");

    %% compute data recycling opinf
    tX_b0_ = tX_b0(1:n_,:);
    [Q,R,P] = qr(tX_b0_,"econ");
    [p,~] = find(P(:,1:n_));  
    [rVn_,~,~] = svd(X_b(:,p),"econ");
    rtO_ = rVn_'*A1*rVn_;

    rtX_b0 = rVn_'*X_b(:,p);
    rtX_b1 = rVn_'*X_b(:,p+1);
    dot_rtX_b = (rtX_b1 - rtX_b0)/dt;
    [rO_,~,~,r_condD] = opinf(dot_rtX_b,rtX_b0,U_b(:,p),is,true);

    r_condsD(j) = r_condD;
   
    r_O_errors(j) = norm(rO_-rtO_,"fro")/norm(rtO_,"fro");
    

    %% compute avg ROM state error
    computeROMStateError = true;
    % computeROMStateError = false
    if computeROMStateError
        Vn_ = Vn(:,1:n_);
        tf = @(tx,u) tO_*tx;
        hO_ = O;
        hf = @(hx,u) hO_*hx;

        sf = @(x,u) sO_*x;
        rf = @(x,u) rO_*x;

        test_type = "train"
        % test_type = "worst-case"

        if test_type == "train"
        x0_r = Vn_'*x0;
        x0_r2 = rVn_'*x0;
        X_t = X_b;
        %% compute avg ROM state error for worst-case initial condition
        elseif test_type == "worst-case"
        [sx0_wc,~,~] = svds(tO_-sO_,1); 
        x0_r = sx0_wc;
        X_t = gen_FOM_data(Vn_*x0_r,u_val,f,nt,N,dt); 
        x0_r2 = rVn_'*Vn_*x0_r;
        else
            error("unknown test_type")
        end
        %%

        t_ROM_state_error(j) = compute_avg_rom_state_error(x0_r,tf,nt,U_b,X_t,Vn_,dt);
        h_ROM_state_error(j) = compute_avg_rom_state_error(x0_r,hf,nt,U_b,X_t,Vn_,dt);
        s_ROM_state_error(j) = compute_avg_rom_state_error(x0_r,sf,nt,U_b,X_t,Vn_,dt);

        r_ROM_state_error(j) = compute_avg_rom_state_error(x0_r2,rf,nt,U_b,X_t,rVn_,dt);  
        
        best_approx_error_POD(j)  = norm(Vn_*Vn_'*X_b - X_b,"fro")/norm(X_b,"fro");
        best_approx_error_mPOD(j) = norm(rVn_*rVn_'*X_b - X_b,"fro")/norm(X_b,"fro");
        
    end    
end

figure
hold on
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"O exact opinf")
semilogy(ns,s_O_errors,'x-', 'LineWidth', 2,'DisplayName',"O standard opinf")
semilogy(ns,r_O_errors,'x-', 'LineWidth', 2,'DisplayName',"O data recycling")
ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")
box on

figure
hold on
semilogy(ns,condsD,'x-', 'LineWidth', 2,'DisplayName',"O exact opinf")
semilogy(ns,s_condsD,'x-', 'LineWidth', 2,'DisplayName',"O standard opinf")
semilogy(ns,r_condsD,'x-', 'LineWidth', 2,'DisplayName',"O data recycling")
ylabel("condition number")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')

legend("show")
box on

% save("data_heat","O_errors","condsD");

if computeROMStateError
    figure
    hold on
    semilogy(ns,h_ROM_state_error,'x-', 'LineWidth', 2,'DisplayName',"exactOpInf", "MarkerSize",10)
    semilogy(ns,t_ROM_state_error,'+:', 'LineWidth', 2,'DisplayName',"intrusive", "MarkerSize",10)
    semilogy(ns,s_ROM_state_error,'o:', 'LineWidth', 2,'DisplayName',"standard OpInf", "MarkerSize",10)
    semilogy(ns,r_ROM_state_error,'+:', 'LineWidth', 2,'DisplayName',"data recycling", "MarkerSize",10)
    semilogy(ns,best_approx_error_POD,'-.', 'LineWidth', 2,'DisplayName',"best approx error POD", "MarkerSize",10)
    semilogy(ns,best_approx_error_mPOD,'-.', 'LineWidth', 2,'DisplayName',"best approx error manipulated POD", "MarkerSize",10)
    ylabel("avg rel error of states","Interpreter","latex", "FontSize",15)
    xlabel("ROM dimension","Interpreter","latex", "FontSize",15)
    set(gca, 'YScale', 'log')
    grid on
    legend("show","Interpreter","latex", "FontSize",12)
    legend("Location","northeast")
    % ylim([1e-17 1e-15])
    box on
    % savefig("figures/rom_state_error_chafee_infante.fig")
    % exportgraphics(gcf,"figures/rom_state_error_chafee_infante.pdf")
end


% %% FOM solver running for one time step
% function x_1 = single_step(x_0,u_0,dt,f)
%     x_1 = x_0 + dt*f(x_0,u_0);
% end

