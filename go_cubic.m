clear all;
close all;

rng(1); % for reproducibility

N = 12;

%% Allen-Cahn equation as in https://doi.org/10.1016/j.cma.2022.115836

A1 = diag(ones(N-1,1),-1) -eye(N);
A1 = (A1+A1')*N^2;
JN3 = power2kron(N,3); % J_N^(3)
A3 = vecwise_kron(eye(N),3)';
A3 = A3*JN3;
B = zeros(N,1);
B(1) = 1;

IN3 = kron2power(N,3); % I_N^(3)
f = @(x,u) A1*x + A3*IN3*vecwise_kron(x,3) + B*u;

x0 = zeros(N,1);
u_val = @(t) 10*(sin(pi*t)+1); % U_val

dt = 1e-5;
t_end = .1;
nt = t_end/dt;

X_b = zeros(N,nt);
U_b = zeros(1,nt);

t = 0;
x = x0;
u = u_val(t);

X_b(:,1) = x0;
U_b(:,1) = u;

for i=2:nt
    x = x + dt*f(x,u);
    t = t + dt;
    u = u_val(t);

    X_b(:,i) = x;
    U_b(:,i) = u;
end

