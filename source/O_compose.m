function O = O_compose(ops,is,qs,k)

% initialize parametric aggregated operator with zeros 
% and add operators ops at the k-th position
%
% example: is = [1 2] , qs = [qA qH], ops = [A_k H_k]
% O = [ A_1 ... A_qA | H_1 ... H_qH , but with all entries zero,
% except for A_k and H_k

% ! no input signal u considered! -> Nu = 0

n = size(ops,1);
% n_is_ = n_is(n,is);
% 
% O_lengths = n_is_.*qs;
% O = zeros(n,sum(O_lengths));
O = O_empty(n,is,qs);


% offset1 = [0 O_lengths(1:end-1)];
% offset2 = (k-1)*n_is_;
% ks = [];
% for j = 1:length(is)
%     % only consider the operators whose affine expansions have at least k
%     % elements
%     if qs(j)>=k 
%         ks = [ks offset1(j)+offset2(j)+(1:n_is_(j))];
%     end
% end

ks = get_ops_indices(n,is,qs,k,is);

O(:,ks) = ops;
