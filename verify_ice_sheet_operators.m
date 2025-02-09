N = 3;
% I = speye(N);
I = eye(N);
D = diag(ones(N-1,1),-1) - I; % not sparse!
K = vecwise_kron(I)'; % not sparse!

%% check third-order term
X3 = rank_suff_basis(N,3);

A3 = D*K*kron(D,K*kron(I,I));
XA3 = A3*vecwise_kron(X3,3);

f3 = @(x) D*(D*x.*x.*x);

N_3 = size(X3,2);
Xf3 = zeros(N,N_3);

for i = 1:N_3
    Xf3(:,i) = f3(X3(:,i));
end

XA3 - Xf3
%%