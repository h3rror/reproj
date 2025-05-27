function [D,A_inds,B_inds] = getOpInfMatrix(tX,U,is)

D = U;
ni = numel(is);
n = size(tX,1);
A_inds = zeros(ni,2);
p = size(U,1);
B_inds = [1 p];

ind_old = p;
for k = 1:ni
    i = is(k);
    % if isa(tX,"int32")
    %     [~,rows] = find(kron2power(n,i));
    %     tXi = vecwise_kron(tX,i);
    %     D = [D; tXi(rows,:)];
    % else
    %     D = [D; kron2power(n,i)*vecwise_kron(tX,i)];
    % end
    D = [D; uniquepower(tX,i)];
    ind = size(D,1);
    A_inds(k,:) = [ind_old+1, ind];
    ind_old = ind;
end
D = sparse(double(D));