function I = kron2power(N,i)
% computes matrix I_N^{(2)} as described in Peherstorfer reprojection
% article

if(nargin < 2)
    i = 2;
end

[~,~,u] = reduced_coordinates(N,i);

N_i = nchoosek(N+i-1,i); % see Peherstorfer's reproj article
I = sparse(1:N_i,u,1);

% test: kk = 3; kron2power(kk)*power2kron(kk)
