function rhs = myRK4(f, deltaT, x)
% Runge Kutta 4th order discretization in time

k1 = deltaT*f(x);
k2 = deltaT*f(x + k1/2);
k3 = deltaT*f(x + k2/2);
k4 = deltaT*f(x + k3);
rhs = x + 1/6*(k1 + 2*k2 + 2*k3 + k4);

end