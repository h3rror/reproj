function kron_h = half_kron(vector)
% kronecker product of a vector with itself without ambiguous entries:
% a_i a_j = a_j a_i
%
% in incremental order: high index entries appear as late as possible

n = numel(vector);

counter = 1;
for i = 1:n
    for j = 1:i
        kron_h(counter) = vector(i)*vector(j);
        counter = counter +1;
    end
end