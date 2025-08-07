clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

% is = [1 2 3];
is = [1 2];
% is = 3;

s = 5;
% n = 16;
% n = 50;
n = 100;
% n = 3;

% qs = s;
qs = s*ones(size(is));

tX = rank_suff_basis(n,is);
ptX = repmat(full(tX),1,1,s);

Theta = magic(s);
cond(Theta)

Thetas{1} = Theta;
Thetas{2} = Theta.^2;
% Thetas{3} = Theta.^3;
% Thetas{2} = Theta;
% Thetas{3} = Theta;

% ns = 1:n;
ns = 10:10:n;
% ns = n;

nn = numel(ns);

conds = zeros(nn,1);
pconds = zeros(nn,1);

qs0 = ones(size(is));

for k = 1:nn
    n_ = ns(k);

    tX_ = O_change(tX,is,qs0,qs0,n_);

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
title("polynomial degree "+is)



