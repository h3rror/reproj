clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

is = 2;

s = 5;
n = 16;

qs = s;

tX = rank_suff_basis(n,is);
ptX = repmat(full(tX),1,1,s);

Theta = magic(s);
cond(Theta)

Thetas{1} = Theta;

ns = 1:n;
% ns = n;

nn = numel(ns);

conds = zeros(nn,1);
pconds = zeros(nn,1);

for k = 1:nn
    n_ = ns(k);

    tX_ = O_change(tX,is,1,1,n_);

    D = getOpInfMatrix(tX_,[],is);
    conds(k) = cond(full(D));

    % ptX_ = O_change(ptX(:,:),is,qs,qs,n_);
    ptX_ = repmat(full(tX_),1,1,s);

    pD = getPOpInfMatrix(ptX_,[],is,Thetas);
    pconds(k) = cond(full(pD));
end


figure
hold on
semilogy(ns,pconds,'x-', 'LineWidth', 2,'DisplayName',"monolithic")
semilogy(ns,conds,'x-', 'LineWidth', 2,'DisplayName',"decoupled")
ylabel("condition number")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")



