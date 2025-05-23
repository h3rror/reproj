function reproj(nList)

if(nargin < 1)
    nList = 1:6;
end

%% generate system matrix
N = 10; % full-model dimension
rng(1); % for reproducibility
D = diag(-logspace(-1, -2, N));
[W, ~] = qr(randn(N, N), 0);
tcA = W'*D*W; % system matrix of time-continuous system
save tcA tcA
dlmwrite('tcA.txt', tcA, 'precision', 16, 'delimiter', ' ');

%% discretize with RK4
deltaT = 1; % time step size
K = 100; % number of time steps
f = @(x)tcA*x;
makeTimeStepRK4 = @(x)myRK4(f, deltaT, x);
% time-discrete matrix A
E = eye(N);
A = E + 1/6*deltaT*(tcA + 2*tcA*(E + deltaT/2*tcA) + 2*tcA*(E + deltaT/2*(tcA*(E + deltaT/2*tcA))) + tcA*(E + deltaT*tcA*(E + deltaT/2*(tcA*(E + deltaT/2*tcA)))));
% A = zeros(N, N);
% for i=1:N
%     A(:, i) = makeTimeStepRK4(E(:, i));
% end

%% set initial condition
% training initial condition is [1, 0, ..., 0]^T
% test initial condition is [1, 1, 0, ..., 0]^T
Xtrain = zeros(N, K);
Xtrain(:, 1) = E(:, 1);
Xtest = zeros(N, K);
Xtest(:, 1) = E(:, 1) + E(:, 2);

%% time step full model
for i=2:K
    Xtrain(:, i) = makeTimeStepRK4(Xtrain(:, i - 1));
    Xtest(:, i) = makeTimeStepRK4(Xtest(:, i - 1));
end

for n=nList
    %% Construct basis and project full-model trajectory
    V = E(:, 1:n);
    XtrainProj = V'*Xtrain;
    XtestProj = V'*Xtest;
    
    %% Re-project
    XtrainReProj = zeros(n, K);
    XtrainReProj(:, 1) = V'*Xtrain(:, 1);
    for i=2:K
        XtrainReProj(:, i) = V'*makeTimeStepRK4(V*XtrainReProj(:, i - 1));
    end
    
    %% Intrusive model reduction
    Ar = V'*A*V;
    
    %% Operator inference without re-proj
    ArOpInf = (XtrainProj(:, 1:K-1)'\XtrainProj(:, 2:K)')';
    norm(Ar - ArOpInf)
    
    %% Operator inference with re-proj
    ArOpInfReProj = (XtrainReProj(:, 1:K-1)'\XtrainReProj(:, 2:K)')';
    norm(Ar - ArOpInfReProj)
    
    %% Time stepping reduced models
    XtestIntMOR = zeros(n, K);
    XtrainIntMOR = zeros(n, K);
    XtestIntMOR(:, 1) = V'*Xtest(:, 1);
    XtrainIntMOR(:, 1) = V'*Xtrain(:, 1);
    XtestOpInf = zeros(n, K);
    XtestOpInf(:, 1) = V'*Xtest(:, 1);
    XtestOpInfReProj = zeros(n, K);
    XtestOpInfReProj(:, 1) = V'*Xtest(:, 1);
    for i=2:K
        XtestIntMOR(:, i) = Ar*XtestIntMOR(:, i - 1);
        XtrainIntMOR(:, i) = Ar*XtrainIntMOR(:, i - 1);
        XtestOpInf(:, i) = ArOpInf*XtestOpInf(:, i - 1);
        XtestOpInfReProj(:, i) = ArOpInfReProj*XtestOpInfReProj(:, i - 1);
    end
    
    %% plot
    marker_inds = 1:10:100;
    figure;
    plot(sqrt(sum(XtestProj.^2, 1)), '--g', 'LineWidth', 2,'MarkerIndices',marker_inds);
    hold on;
    plot(sqrt(sum(XtestIntMOR.^2, 1)), '-ok', 'LineWidth', 2,'MarkerIndices',marker_inds);
    hold on;
    plot(sqrt(sum(XtestOpInf.^2, 1)), '-sr', 'LineWidth', 2,'MarkerIndices',marker_inds);
    hold on;
    plot(sqrt(sum(XtestOpInfReProj.^2, 1)), '-mx', 'LineWidth', 2,'MarkerIndices',marker_inds);
    xlabel('time step k');
    ylabel('2-norm of states');
    legend('projected', 'intrusive model reduction', 'OpInf, w/out re-proj', 'OpInf, re-proj');
    title(['Dimension n = ', num2str(n)]);
    axis([-Inf Inf 0 1.6]);
    
end

end

function rhs = myRK4(f, deltaT, x)
% Runge Kutta 4th order discretization in time

k1 = deltaT*f(x);
k2 = deltaT*f(x + k1/2);
k3 = deltaT*f(x + k2/2);
k4 = deltaT*f(x + k3);
rhs = x + 1/6*(k1 + 2*k2 + 2*k3 + k4);

end