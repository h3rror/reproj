clear all;
close all;

rng(1); % for reproducibility

% synthetic example (11) in
% https://www.sciencedirect.com/science/article/pii/S0045782522007927?casa_token=AXNbPCfGtzoAAAAA:6IFh6fogD6SLiCtQxXdhGaCwHohuZCcKVcBsBB1h7Y55JjrUkMinEEL5mtDtx2jTEOaqj9IOoj8

N = 128;
% N = 10;
N_2 = N*(N+1)/2;

B = rand(N,1);
F = rand(N,N_2); % problem: leads to unstable simulations
%% better F: energy-conserving
F_s = rand(N,N,N);
F_s2 = F_s - permute(F_s,[2 1 3]);
F = reduced_convection_operator(F_s2(:,:));
%% FOM snapshot generation for ROM basis construction
A_s = rand(N,N);
mus = (1:10)/10;
% mus = 1;

T = 1;
dt = 1e-3;
nT = T/dt;
nmu = numel(mus);
X_b = zeros(N,nT,nmu);

for i = nmu
    mu = mus(i);
    A = - mu*(A_s+A_s'+ 2*N*eye(N));

    U_1 = 2*rand(1,nT);
    x_0 = rand(N,1);
    X_b(:,1,i) = x_0;

    x = x_0;
    for j = 1:nT
        x_2 = vectorwise_halfkron(x);
        u = U_1(:,j);
        x = x + dt*(A*x + F*x_2 + B*u);
        X_b(:,j+1,i) = x;
    end
end

%% construct POD basis
[Phi,S,~] = svd(X_b(:,:),'econ');
% Phi = eye(10);

ns = [2 4 6 8 10];
% ns = 10;
nn = numel(ns);

e_train = zeros(nn,1);
h_e_train = zeros(nn,1);

A_error = zeros(nn,1);
F_error = zeros(nn,1);
B_error = zeros(nn,1);

Mt = 1;
U_train = U_1;
x_0train = x_0;

for k = 1:nn
    n = ns(k);
    V = Phi(:,1:n);
    tB = V'*B;
    tF = V'*F*kron2power(N)*kron(V,V)*power2kron(n);
    mu = .7;
    A = - mu*(A_s+A_s'+ 2*N*eye(N));
    tA = V'*A*V;

    tX_train = zeros(n,nT,Mt);
    X_train = zeros(N,nT,Mt);

    hX_train = zeros(n,nT,Mt);

    %% intrusive ROM training error
    for i = 1:Mt
        U = U_train(:,:,i);
        x_0 = x_0train(:,i);
        tx_0 = V'*x_0;

        X_train(:,1,i) = x_0;
        tX_train(:,1,i) = tx_0;

        x = x_0;
        tx = tx_0;
        for j = 1:nT
            x_2 = vectorwise_halfkron(x);
            u = U(:,j);
            x = x + dt*(A*x + F*x_2 + B*u);
            X_train(:,j+1,i) = x;

            tx_2 = vectorwise_halfkron(tx);
            u = U(:,j);
            %% naive ROM
            xr = V*tx;
            xr_2 = vectorwise_halfkron(xr);
            tx2 = V'*(xr + dt*(A*xr + F*xr_2 + B*u));
            %%
            tx = tx + dt*(tA*tx + tF*tx_2 + tB*u);
            tX_train(:,j+1,i) = tx;
        end

        error = norm(V*tX_train(:,:,i) - X_train(:,:,i),"fro")/norm(X_train(:,:,i),"fro");
        e_train(k) = e_train(k) + error;
    end
    %% operator inference
    bX = V'*X_train(:,:);
    bX0 = bX(:,1:end-1);
    bX1 = bX(:,2:end);
    bX0_2 = vectorwise_halfkron(bX0);
    D = [bX0; bX0_2; U_train(:,:)];
    bX_dot = (bX1-bX0)/dt;          % watch out: don't take time derivatives across trajectories!
    
    O = (D'\bX_dot')';
    hA = O(:,1:n);
    n_2 = n*(n+1)/2;
    hF = O(:,n+1:n+n_2);
    hB = O(:,n+n_2+1:end);

    A_error(k) = norm(tA-hA);
    F_error(k) = norm(tF-hF);
    B_error(k) = norm(tB-hB);

    for i = 1:Mt
        U = U_train(:,:,i);
        x_0 = x_0train(:,i);
        hx_0 = V'*x_0;

        X_train(:,1,i) = x_0;
        hX_train(:,1,i) = hx_0;

        x = x_0;
        hx = tx_0;
        for j = 1:nT
            % x_2 = vectorwise_halfkron(x);
            % u = U(:,j);
            % x = x + dt*(A*x + F*x_2 + B*u);
            % X_train(:,j+1,i) = x;

            hx_2 = vectorwise_halfkron(hx);
            u = U(:,j);
            %% naive ROM
            % xr = V*tx;
            % xr_2 = vectorwise_halfkron(xr);
            % tx2 = V'*(xr + dt*(A*xr + F*xr_2 + B*u));
            %%
            hx = hx + dt*(hA*hx + hF*hx_2 + hB*u);
            hX_train(:,j+1,i) = hx;
        end

        error = norm(V*hX_train(:,:,i) - X_train(:,:,i),"fro")/norm(X_train(:,:,i),"fro");
        h_e_train(k) = h_e_train(k) + error;
    end
end

figure
semilogy(ns,e_train, "ks--","DisplayName","intrusive")
hold on
semilogy(ns,h_e_train, "ob","Displayname", "OpInf w/o reg")
ylim([1e-5 1])

figure
semilogy(ns,A_error, "rs--","DisplayName","A")
hold on
semilogy(ns,F_error,"bx-","DisplayName","F")
semilogy(ns,B_error,"g+:","DisplayName","B")
ylabel("operator error")
xlabel("ROM dimension")

