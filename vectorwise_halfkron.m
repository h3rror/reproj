function A_kron = vectorwise_halfkron(A)
% vectorwise kronecker product after removing ambiguous entries: 
% e.g.,  a_i a_j = a_j a_i

r = size(A,1);
K = size(A,2);

A_kron = zeros(r*(r+1)/2,K);
for k = 1:K
    a_j = A(:,k);
    A_kron(:,k) = half_kron(a_j);
end