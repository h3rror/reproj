clear all;
close all;

rng(1); % for reproducibility

% synthetic example (11) in
% https://www.sciencedirect.com/science/article/pii/S0045782522007927?casa_token=AXNbPCfGtzoAAAAA:6IFh6fogD6SLiCtQxXdhGaCwHohuZCcKVcBsBB1h7Y55JjrUkMinEEL5mtDtx2jTEOaqj9IOoj8

N = 128;
N_2 = N*(N+1)/2;

B = rand(N,1);
F = rand(N,N_2); % problem: leads to unstable simulations
%% better F: energy-conserving
F_s = rand(N,N,N);
F_s2 = F_s - permute(F_s,[2 1 3]);
F = reduced_convection_operator(F_s2(:,:));
%%
A_s = rand(N,N);
mus = (1:10)/10;
% mus = 1;

T = 1;
dt = 1e-3;
nT = T/dt;
nmu = numel(mus);
X_train = zeros(N,nT,nmu);

for i = nmu
    mu = mus(i);
    A = - mu*(A_s+A_s'+ 2*N*eye(N));

    U_1 = 2*rand(1,nT);
    x_0 = rand(N,1);
    X_train(:,1,i) = x_0;

    x = x_0;
    for j = 1:nT
        x_2 = vectorwise_halfkron(x);
        u = U_1(:,j);
        x = x + dt*(A*x + F*x_2 + B*u);
        X_train(:,j+1,i) = x;
    end
end

%% construct POD basis
[U,S,~] = svd(X_train(:,:),'econ');

ns = [2 4 6 8 10];
nn = numel(ns);

e_train = zeros(nn,1);

for k = nn
    n = ns(k);
    V = U(:,1:n);
    tB = V'*B;
    tF = V'*F*vectorwise_halfkron(V);

    tX_train = zeros(n,nT,nmu);

    %% intrusive ROM training error
    for i = nmu
        mu = mus(i);
        A = - mu*(A_s+A_s'+ 2*N*eye(N));
        tA = V'*A*V;

        U_1 = 2*rand(1,nT);
        x_0 = rand(N,1);
        tx_0 = V'*x_0;
        tX_train(:,1,i) = x_0;

        x = x_0;
        for j = 1:nT
            x_2 = vectorwise_halfkron(x);
            u = U_1(:,j);
            x = x + dt*(A*x + F*x_2 + B*u);
            tX_train(:,j+1,i) = x;
        end
    end

    
    %%
end

