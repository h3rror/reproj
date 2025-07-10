clear all;
close all;

rng(1); % for reproducibility

addpath('source/');

%%
n = 3;
is = [2];

% basis0 = rank_suff_basis(n,is);
basis0 = full(rank_suff_basis(n,is));
basis1 = basis0;
K = size(basis1,2);

%% options
b_rand = true;
% b_rand = false;
b_trunc = true;

%% user-specified tolerances
t_sigma = 1e-7;
t_basis = 1e-3;

%% random weighting of snapshot data
if b_rand
    R = rand(K,K);
    [Q,~,~] = svd(R);
    basis1 = basis1*Q;
end
%% truncate basis
if b_trunc
    n1 = n; % <= K
    % n1 = K-1; % <= K
    basis1 = basis1(:,1:n1); 
end
%% generate snapshot data
P0 = uniquepowers(basis0,is);
P1 = uniquepowers(basis1,is);

A = randi(n^2,n,K);
f = @(X) A*uniquepowers(X,is);
Pdot1 = f(basis1);

%% detect rank deficiency and fix it
if b_trunc
    [U,S,~] = svd(P1);
    diag(S)
    n11 = sum(diag(S)>t_sigma)
    S_trunc = S(1:n11,1:n11);
    P11 = U(:,1:n11);
    Pdot11 = Pdot1*P1'*P11*inv(S_trunc)^2;
    recovery = vecwise_2norm(P11*P11'*P0)./vecwise_2norm(P0);
    sorted = sort(recovery);
    treshold_ = max(t_basis,sorted(K-n11))
    % basis2 = basis0(:,vecwise_2norm(P0 - P11*P11'*P0)>t_basis);
    basis2 = basis0(:,recovery <= treshold_);
    P2 = uniquepowers(basis2,is);

    Pdot2 = f(basis2);
    Pdot = [Pdot11 Pdot2];
    P = [P11 P2];
end
%% exact OpInf

Ahat = Pdot/P
norm(Ahat-A)




