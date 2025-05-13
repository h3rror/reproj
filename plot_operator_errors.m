clear all;
close all;

figure(1) % operator errors
hold on
figure(2) % condition numbers
hold on


load("data_chafee_infante","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"Chafee-Infante")
figure(2)
semilogy(ns,condsD,'x-', 'LineWidth', 2,'DisplayName',"Chafee-Infante")

load("data_icesheet","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"ice sheet")
figure(2)
semilogy(ns,condsD,'x-', 'LineWidth', 2,'DisplayName',"ice sheet")

load("data_burgers","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'x-', 'LineWidth', 2,'DisplayName',"Burgers'")
figure(2)
semilogy(ns,condsD,'x-', 'LineWidth', 2,'DisplayName',"Burgers'")

figure(1)
ylabel("relative operator error")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")

figure(2)
ylabel("condition number")
xlabel("ROM dimension")
set(gca, 'YScale', 'log')
grid on
legend("show")

