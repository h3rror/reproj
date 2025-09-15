function X = simulate(x0,dt,nt,single_step)

N = size(x0,1);
X = zeros(N,nt+1);

t = 0;
x = x0;

X(:,1) = x0;

for i=1:nt
    x = single_step(x);
    t = t + dt;

    X(:,i+1) = x;
end