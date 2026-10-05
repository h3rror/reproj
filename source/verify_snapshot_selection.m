

n = n_;

rel_proj_errors = vecwise_2norm(Vn_*Vn_'*X_b-X_b)./vecwise_2norm(X_b);

figure; plot(rel_proj_errors,"DisplayName", "all snapshots - full set POD")

hold on
% for jj=1:n
for jj=n:n
    % for jj=31:33
    % plot(rel_proj_errors,"s","MarkerIndices",selects{jj},"DisplayName","recycled snapshots "+num2str(jj))
    plot(rel_proj_errors,"s","MarkerIndices",p,"DisplayName","recycled snapshots "+num2str(jj))
end

for jj = 1:n
    [rVn_,~,~] = svd(X_b(:,p(1:jj)),"econ");
    rel_proj_errors_j = vecwise_2norm(rVn_*rVn_'*X_b-X_b)./vecwise_2norm(X_b);

    plot(rel_proj_errors_j,"DisplayName", "all snapshots - slect POD - iteration "+num2str(jj))
end
