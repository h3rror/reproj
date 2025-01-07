close all
clear all
% plot norm of inverse of phi(B) and condition number of phi(B) for several
% ROM dimensions r

% rs = 1:12;
% rs = 100;
% rs = 40;
rs = 1:100;

% no_rs = numel(rs);
r_max = max(rs);
inv_norms = zeros(r_max,1);
conds = zeros(r_max,1);

for r = rs
    B = basis_lin_quad(r);
    A = [B; vectorwise_halfkron(B)];
    inv_A = inv(A);
    inv_norms(r) = norm(inv_A);
    conds(r) = cond(A);
end

figure
plot(rs,conds,'x:')
hold on
plot(rs,inv_norms,'x:')

legend("cond$(\varphi(B))$","$||\varphi(B)^{-1}||$",'Location','northwest','Interpreter','latex')

plot(rs,rs, 'DisplayName',"r")
plot(rs,sqrt(rs),"displayname","$\sqrt(r)$")

set(gca,'Yscale','log');
set(gca,'Xscale','log');

ylabel("value")
xlabel("ROM dimension r")
