function dt = dt_estimate(X_b,U_b,Vn1,dt)

Vn1 = Vn(:,1);
X1 = X_b(:,1:end-1);
X2 = X_b(:,2:end);
dot_X = (X2-X1)/dt;
U_b1 = U_b(:,1:end-1);

vecwise_norm = @(M) sqrt(sum(M.^2,1));

dt = 1/max(vecwise_norm(Vn1'*dot_X)./vecwise_norm(getOpInfMatrix(Vn1'*X1,U_b1,is)));
