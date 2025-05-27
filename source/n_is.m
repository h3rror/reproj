function n_is = n_is(ns,is)
% compute n_i as defined in https://epubs.siam.org/doi/epdf/10.1137/19M1292448
% for array of is
% and arrays of ns

n_i = @(n,i) nchoosek(n+i-1,i);

% n_is = [];
% for i = is
%     n_is = [n_is n_i(n,i)];
% end

nn = numel(ns);
ni = numel(is);
n_is = zeros(nn,ni);

for j=1:nn
    n = ns(j);
    for k =1:ni
        i = is(k);
        n_is(j,k) = n_i(n,i);
    end
end

