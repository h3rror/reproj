function I = kron2power(N)
% computes matrix I_N^{(2)} as described in Peherstorfer reprojection
% article

% I = eye(N^2);
[~,~,u] = reduced_coordinates(N);
% I = I(u,:);

% test: kk = 3; kron2power(kk)*power2kron(kk)


% sparse
N_2 = N*(N+1)/2;
I = sparse(1:N_2,u,1);