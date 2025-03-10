clear all;
close all;

rng(1); % for reproducibility

figure
hold on

%% burgers

% title("Burgers")
% is = [1 2];
% % is_c = 0:max(is);
% 
% n_trial = 5*10^4;
% ns_trial = 1:10;

%% chafee-infantee
% title("Chafee-Infante")
% is = [1 3];
% 
% n_trial = 10^7;
% ns_trial = 1:12;
% ylim([1 10^8])

%% reaction-diffusion
title("Reaction-diffusion")

is = [1 2 3];
n_trial = 5*10^3;
ns_trial = 1:10;


ns = 1:max(ns_trial)+10;

n_is_ = n_is(ns,is);
n_snaps = sum(n_is_,2);


plot(ns,n_snaps,'x-','LineWidth', 2,'DisplayName',"rank-suff")
plot(ns_trial,n_trial*ones(size(ns_trial)),'x-','LineWidth', 2,'DisplayName',"trial and error")

set(gca, 'YScale', 'log')
ylabel("number of snapshots")
xlabel("ROM dimension")

% title("training error")
grid on
legend("show")
