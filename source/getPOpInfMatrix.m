function [D,A_inds,B_inds] = getPOpInfMatrix(tX,U,is,Thetas)

% warning("not handling inputs u yet")

% D = U;
D = [];
ni = numel(is);
n = size(tX,1);
A_inds = zeros(ni,2);
p = size(U,1);
B_inds = [1 p];

s = size(tX,3);

ind_old = p;
for k = 1:ni
    i = is(k);
    Theta = Thetas{k};
    D_column = [];
    for j =1:s
        theta = Theta(j,:);
        D_column = [D_column; kron(theta,uniquepower(tX(:,:,j),i)')];
    end
    D = [D, D_column];
    ind = size(D,2);
    A_inds(k,:) = [ind_old+1, ind];
    ind_old = ind;
end

D = D';