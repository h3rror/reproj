N = 8;
% I = speye(N);
I = eye(N);
D = diag(ones(N-1,1),-1) - I; % not sparse!
K = vecwise_kron(I)'; % not sparse!

%% check third-order term
% X3 = rank_suff_basis(N,3);
% 
% A3 = D*K*kron(D,K*kron(I,I));
% XA3 = A3*vecwise_kron(X3,3);
% 
% f3 = @(x) D*(D*x.*x.*x);
% 
% N_3 = size(X3,2);
% Xf3 = zeros(N,N_3);
% 
% for i = 1:N_3
%     Xf3(:,i) = f3(X3(:,i));
% end
% 
% XA3 - Xf3
%% check eighth-order term
X8 = rank_suff_basis(N,8);

A8 = D*K*kron(D,K*kron(D,K*kron(D,K*kron(I,K*kron(I,K*kron(I,K*kron(I,I)))))));
XA8 = A8*vecwise_kron(X8,8);

f8 = @(x) D*(x.^5.*(D*x).^3);

N_8 = size(X8,2);
Xf8 = zeros(N,N_8);

for i =1:N_8
    Xf8(:,i) = f8(X8(:,i));
end

XA8 - Xf8
