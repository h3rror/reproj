function O_k = extract(O,is,qs,k,is_select)

% extract k-th suboperator from parametric aggregated operator
%
% example: is = [1 2] , qs = [qA qH], is_select = [2]
% so O = [ A_1 ... A_qA | H_1 ... H_qH , then output [H_k]