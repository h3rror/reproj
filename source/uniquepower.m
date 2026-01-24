function Xi = uniquepower(X,i,u)
% computes vector x^i for all columns x of X without memory-expensive
% detour via kronecker products of x

[N,K] = size(X);

if(nargin < 3)
    [s,b,u] = reduced_coordinates(N,i);
end

Ni = numel(u);
% subs_s = zeros(Ni);
Xi = zeros(Ni,K);

for j = 1:Ni
    subs_s = ind2subs(N*ones(1,i),u(j));
    Xi(j,:) = prod(X(subs_s,:),1);
end


