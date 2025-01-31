function J = power2kron(N,i)
% computes matrix J_N^{(2)} as described in Peherstorfer reprojection
% article

if(nargin < 2)
    i = 2;
end

% N_2 = N*(N+1)/2;
% 
% J = eye(N_2);
s = reduced_coordinates(N,i);
% J = J(s,:);

% test: kk = 3; kron2power(kk)*power2kron(kk)

% sparse
J = sparse(1:N^i,s,1);
