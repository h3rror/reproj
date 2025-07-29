function [O,A_inds,B_inds,condD] = p_opinf(dot_tX,tX,U,is,Thetas,verbose)
    % dot_tX: r x K x s    ROM state time derivative snapshots per parameter
    % tX: r x K x s    ROM state snapshots per parameter
    % is: array of polynomial degrees
    % Thetas: cell array of numel(is) matrices Theta_i, 
    % each of dimension s x q_{A_i}

    % O: aggregated ROM operator
    % A_inds: length(is) x 2 array of indices such that
    %         A_is(k) = O(:,A_inds(k,:)) for all k = 1,...,length(is)
    % B_inds: 1 x 2 array of indices such that B = O(:,B_inds);

    % D = [];

    
if(nargin < 6)
    verbose = false;
end

    [D,A_inds,B_inds] = getPOpInfMatrix(tX,U,is,Thetas);

    if verbose==true
        condD = cond(full(D));
    else
        condD = false;
    end

   O = (D'\dot_tX(:,:)')';
end