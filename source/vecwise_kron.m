function A_kron = vecwise_kron(A,i)

if(nargin < 2)
    i = 2;
end

[r,K] = size(A);

if isa(A,"int32")
    A_kron = zeros(r^i,K, "int32");
else
    A_kron = zeros(r^i,K);
end

for j=1:K
    a_j = A(:,j);
    a_kron = a_j;
    for l=1:i-1
        a_kron = kron(a_kron,a_j);
    end
    A_kron(:,j) = a_kron;
end

if issparse(A)
    A_kron = sparse(A_kron);
end