clear all;
close all;

rng(1); % for reproducibility

figure
hold on

cases = ["heat", ...
    ... "Burgers'", ...
  "Burgers'", "Chafee-Infante", "Diffusion-reaction", "FitzHugh-Nagumo", "Euler", ...
  "periodic Burgers'", "Kuramoto-Sivashinsky"];
cases = cases(end:-1:1)

n_snapshots = [10^3 9; ....
    ....10^5 120; 
    5*10^4 136; 10^7 455; 5*10^3 286; 3600 135; 64*10^3 528; ...
    7500 65; 2.7e4 324];
n_snapshots = n_snapshots(end:-1:1,[2 1])
b = barh(cases,n_snapshots);

xtips1 = b(1).YEndPoints + 0.3;
xtips1 = 1.5*b(1).YEndPoints;
ytips1 = b(1).XEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'VerticalAlignment','middle')
% text(xtips1,ytips1,labels1,'VerticalAlignment','right')

% xtips2 = b(2).YEndPoints + 0.3;
xtips2 = 1.5*b(2).YEndPoints;
ytips2 = b(2).XEndPoints;
% labels2 = string(b(2).YData);
labels2 = compose("%1.1e",b(2).YData);
text(xtips2,ytips2,labels2,'VerticalAlignment','middle')

set(gca, 'XScale', 'log')
% grid on
xlabel("number of snapshots")
xlim([1. 2e8])

% legend({'trial and error', 'rank-ensuring'})
legend({'rank-ensuring', 'trial and error'})