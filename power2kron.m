function J = power2kron(N)
% computes matrix J_N^{(2)} as described in Peherstorfer reprojection
% article

N_2 = N*(N+1)/2;

J = eye(N_2);
s = reduced_coordinates(N);
J = J(s,:);

% test: kk = 3; kron2power(kk)*power2kron(kk)