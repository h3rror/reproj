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
% b_rand = true;
b_rand = false;
b_trunc = true;

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
    [P_ortho,S,~] = svd(P1);
    diag(S)
    % P2 = P_ortho(:,n1+1:end);
    P2 = P0(:,n1+1:end); % cheat
    basis2 = basis0*inv(P0)*P2;

    Pdot2 = f(basis2);
    Pdot = [Pdot1 Pdot2];
end
%% exact OpInf

Ahat = Pdot/P0
norm(Ahat-A)




