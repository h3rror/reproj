%% rank sufficient operator inference
A_error_suff = zeros(nn,1);
F_error_suff = zeros(nn,1);
B_error_suff = zeros(nn,1);

for i = 1:nn
    n = ns(i);
    tX_0_suff = [basis_lin_quad(n), zeros(n,1)]; % zeros to account for input U
    n_2 = n*(n+1)/2;
    U_suff = [zeros(1,n+n_2) 1];

    V = Phi(:,1:n);
    X_0_suff = V*tX_0_suff;
    X_1_suff = zeros(N,n+n_2,1);

    A = - mu*(A_s+A_s'+ 2*N*eye(N));

    for j = 1:n+n_2+1
        x = X_0_suff(:,j);
        u = U_suff(:,j);
        X_1_suff(:,j) = x + dt*(A*x + F*kron2power(N)*kron(x,x) + B*u);
    end

    tX_0_suff = V'*X_0_suff;
    tX_1_suff = V'*X_1_suff;

    tX_dot_suff = (tX_1_suff - tX_0_suff)/dt;

    tX_0_suff_2 = kron2power(n)*vectorwise_kron(tX_0_suff);
    D_suff = [tX_0_suff; tX_0_suff_2; U_suff];
    
    O = (D_suff'\tX_dot_suff')';
    hA = O(:,1:n);
    n_2 = n*(n+1)/2;
    hF = O(:,n+1:n+n_2);
    hB = O(:,n+n_2+1:end);

    tB = V'*B;
    tF = V'*F*kron2power(N)*kron(V,V)*power2kron(n);
    tA = V'*A*V;

    A_error_suff(i) = norm(tA-hA);
    F_error_suff(i) = norm(tF-hF);
    B_error_suff(i) = norm(tB-hB);
end

figure(4)
semilogy(ns,A_error_suff, "rs--","DisplayName","A")
hold on
semilogy(ns,F_error_suff,"bx-","DisplayName","F")
semilogy(ns,B_error_suff,"g+:","DisplayName","B")
ylabel("operator error")
xlabel("ROM dimension")
legend("show")