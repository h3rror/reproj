function O = compose(ops,is,qs,k)

% initialize parametric aggregated operator with zeros 
% and add operators ops at the k-th position
%
% example: is = [1 2] , qs = [qA qH]
% O = [ A_1 ... A_qA | H_1 ... H_qH , but with all entries zero,
% except for A_k and H_k