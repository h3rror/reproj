function O_new = O_change(O,is,qs,qs_new,n_new)
% transform aggregated operator with qs affine terms into 
% operator with qs_new affine terms 
% and dimension n_new

if (nargin<5)
    n = size(O,1);
    n_new = n;
end

% n = size(O,1);
O_new = O_empty(n_new,is,qs_new);
for k=1:max(qs)
    is_select = is(qs_new>=k);
    ops = O_extract(O,is,qs,k,is_select,n_new);
    O_new = O_new + O_compose(ops,is,qs_new,k);
end