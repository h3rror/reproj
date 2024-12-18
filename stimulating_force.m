function f = stimulating_force(V,t,dt,B,phi)
% V: current velocity/state
% t: current time
% dt: time step size
% B: rank sufficient ROM coefficient basis
% phi: ROM basis

r = size(B,1);
B = [B zeros(r,1)]; % add a zero column for last time step r^*+1


ind = round(t/dt);

f = 1/dt*(phi*B(:,ind+2) - V);


