function [rtX_b0,rVn_,p] = data_recycling(tX_b0_, X_b, Vn_)
% select snapshots for data recycling without projection error tolerance

[n_,K] = size(tX_b0_);

% [Q,R,P] = qr(tX_b0_,"econ");
[Q,R,P] = qr(X_b,"econ");  % minimizes largest snapshot projection error
[p,~] = find(P(:,1:n_));
[rVn_,~,~] = svd(X_b(:,p),"econ");

%% verify visually minimization of largest snapshot projection error
verify_snapshot_selection

rtX_b0 = rVn_'*X_b(:,p);
