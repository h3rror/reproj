clear all;

N = 6;
% I = speye(N);
% D = diag(ones(N-1,1),-1) - I; % not sparse!
% D = sparse(D);
% K = vecwise_kron(I)'; 

%% check eighth-order term
X8 = rank_suff_basis(N,8);

% A8 = D*K*kron(D,K*kron(D,K*kron(D,K*kron(I,K*kron(I,K*kron(I,K*kron(I,I)))))));
% XA8 = A8*vecwise_kron(X8(:,1:600),8);

X8_kron_1 = vecwise_kron(X8(:,1:600),8);
X8_kron_2 = vecwise_kron(X8(:,601:end),8);
X8_kron = [X8_kron_1 X8_kron_2];
