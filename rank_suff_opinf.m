%% rank sufficient operator inference
A_error_suff = zeros(nn,1);
F_error_suff = zeros(nn,1);
B_error_suff = zeros(nn,1);

% for i = 1:nn
for i = nn % only run once for largest n
    n = ns(i);
    tX_0_suff = [zeros(n,1), basis_lin_quad(n)]; % zeros to account for input U
    n_2 = n*(n+1)/2;
    U_suff = [1, zeros(1,n+n_2)];

    V = Phi(:,1:n);
    X_0_suff = V*tX_0_suff;
    X_1_suff = zeros(N,n+n_2,1);

    A = - mu*(A_s+A_s'+ 2*N*eye(N));

    I_Npow2 = kron2power(N); % peherstorfer reproj notation

    for j = 1:n+n_2+1
        x = X_0_suff(:,j);
        u = U_suff(:,j);
        X_1_suff(:,j) = x + dt*(A*x + F*I_Npow2*kron(x,x) + B*u);
    end

    tX_0_suff = V'*X_0_suff;
    tX_1_suff = V'*X_1_suff;

    tX_dot_suff = (tX_1_suff - tX_0_suff)/dt;

    tX_0_suff_2 = kron2power(n)*vectorwise_kron(tX_0_suff); % should work to do this only once thanks to incrementality-preserving ordering

    %% compute intrusive operators (has actually already been done before)
    tB = V'*B;
    tF = V'*F*I_Npow2*kron(V,V)*power2kron(n);
    tA = V'*A*V;
end

for i = 1:nn
    n = ns(i);
    n_star = 1+n+n*(n+1)/2;

    tX_0_suff_n = tX_0_suff(1:n,1:n_star);
    tX_1_suff_n = tX_1_suff(1:n,1:n_star);

    tX_dot_suff_n = tX_dot_suff(1:n,1:n_star);

    n_2 = n*(n+1)/2;
    tX_0_suff_2_n = tX_0_suff_2(1:n_2,1:n_star); % maybe doesn't work after all
    % tX_0_suff_2_n = kron2power(n)*vectorwise_kron(tX_0_suff_n);

    D_suff_n = [tX_0_suff_n; tX_0_suff_2_n; U_suff(1:n_star)];
    
    O = (D_suff_n'\tX_dot_suff_n')';
    hA = O(:,1:n);
    n_2 = n*(n+1)/2;
    hF = O(:,n+1:n+n_2);
    hB = O(:,n+n_2+1:end);

    tB_n = tB(1:n,:);
    tA_n = tA(1:n,1:n);
    tF_n = tF(1:n,1:n_2); % only works thanks to incrementality-preserving kron2power/power2kron operators

    A_error_suff(i) = norm(tA_n-hA);
    F_error_suff(i) = norm(tF_n-hF);
    B_error_suff(i) = norm(tB_n-hB);
end

figure
semilogy(ns,A_error_suff, "rs-","DisplayName","A")
hold on
semilogy(ns,F_error_suff,"bx-","DisplayName","F")
semilogy(ns,B_error_suff,"g+-","DisplayName","B")
ylabel("operator error")
xlabel("ROM dimension")
legend("show")