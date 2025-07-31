function ops = O_extract(O,is,qs,k,is_select,n_select)

% extract k-th suboperator from parametric aggregated operator
%
% example: is = [1 2] , qs = [qA qH],  so O = [ A_1 ... A_qA | H_1 ... H_qH ]
% is_select = [1 2] --> output [A_k H_k]
% is_select = [2]   --> output [H_k]
%
% ! no input signal u considered! -> Nu = 0
% is_select gets sorted in increasing order internally, so [2 1] has the
% same output as [1 2]

if (nargin < 6)
    n_select = n;
end

n = size(O,1);
ks = get_ops_indices(n,is,qs,k,is_select,n_select);

ops = O(1:n_select,ks);