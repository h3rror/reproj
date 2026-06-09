function [t_ROM_state_error_j,tX] = compute_avg_rom_state_error(tx0,tf,nt,U_b,X_b,Vn,dt)

n = size(tx0,1);

t = 0;
tx = tx0;
u = U_b(:,1);

tX = zeros(n,nt+1);
tX(:,1) = tx;

for i=1:nt
    tx = single_step(tx,u,dt,tf);

    t = t + dt;
    u = U_b(:,i);

    tX(:,i+1) = tx;
end

t_ROM_state_error_j = norm(Vn*tX - X_b,"fro")/norm(X_b,"fro");