function norms_vec = vecwise_2norm(M)

% compute 2-norm for all columns of M

norms_vec = sqrt(sum(M.^2,1));