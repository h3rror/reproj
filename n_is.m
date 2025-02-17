function n_is = n_is(n,is)
% compute n_i as defined in https://epubs.siam.org/doi/epdf/10.1137/19M1292448
% for array of is

n_i = @(n,i) nchoosek(n+i-1,i);

n_is = [];
for i = is
    n_is = [n_is n_i(n,i)];
end
