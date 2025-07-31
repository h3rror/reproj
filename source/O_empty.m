function O = O_empty(n,is,qs)

% construct all-zero operator of correct dimension
n_is_ = n_is(n,is);

O_lengths = n_is_.*qs;
O = zeros(n,sum(O_lengths));