function [X_b,U_b] = gen_FOM_data(x0,u_val,f,nt,N,dt)

X_b = zeros(N,nt+1);
U_b = zeros(0,nt+1);

t = 0;
x = x0;
u = u_val(t);

X_b(:,1) = x0;
U_b(:,1) = u;

for i=1:nt
    x = x + dt*f(x,u);
    t = t + dt;
    u = u_val(t);

    X_b(:,i+1) = x;
    U_b(:,i+1) = u;
end