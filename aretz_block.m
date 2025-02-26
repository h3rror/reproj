function Q_ij = aretz_block(i,j,n,l)
% compute block matrix of permuted P as suggested by Aretz

xs = nnz_xs(i,n,l);
alphas = nnz_xs(j,n,l);

n_xs = size(xs,2);
n_alphas = size(alphas,2);

Q_ij = zeros(n_xs,n_alphas);

for k = 1:n_alphas
    alpha = alphas(:,k);
    Q_ij(:,k) = prod(xs.^alpha,1);
end