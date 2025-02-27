function [O,A_inds,B_inds] = opinf(dot_tX,tX,U,is)
    % dot_tX: ROM state time derivative snapshots
    % tX: ROM state snapshots
    % is: array of polynomial degrees

    % O: aggregated ROM operator
    % A_inds: length(is) x 2 array of indices such that
    %         A_is(k) = O(:,A_inds(k,:)) for all k = 1,...,length(is)
    % B_inds: 1 x 2 array of indices such that B = O(:,B_inds);

    % D = [];
    D = U;
    ni = numel(is);
    n = size(tX,1);
    A_inds = zeros(ni,2);
    p = size(U,1);
    B_inds = [1 p];

    ind_old = p;
    for k = 1:ni
        i = is(k);
        if isa(tX,"int32")
            [~,rows] = find(kron2power(n,i));
            tXi = vecwise_kron(tX,i);
            D = [D; tXi(rows,:)]; 
        else
            D = [D; kron2power(n,i)*vecwise_kron(tX,i)];
        end
        ind = size(D,1);
        A_inds(k,:) = [ind_old+1, ind];
        ind_old = ind;
    end
    D = sparse(double(D));

   O = (D'\dot_tX')';
end