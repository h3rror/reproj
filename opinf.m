function [O,A_inds,B_inds] = opinf(dot_tX,tX,U,is)
    % dot_tX: ROM state time derivative snapshots
    % tX: ROM state snapshots
    % is: array of polynomial degrees

    % O: aggregated ROM operator
    % A_inds: length(is) x 2 array of indices such that
    %         A_is(k) = O(:,A_inds(k,:)) for all k = 1,...,length(is)
    % B_inds: 1 x 2 array of indices such that B = O(:,B_inds);

    % D = [];

    [D,A_inds,B_inds] = getOpInfMatrix(tX,U,is);

   O = (D'\dot_tX')';
end