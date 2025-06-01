clear all;
close all;

figure(1) % operator errors
hold on
figure(2) % condition numbers
hold on


load("data/data_chafee_infante","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'^-', 'LineWidth', 2,'DisplayName',"Chafee-Infante", "MarkerSize",10)
figure(2)
semilogy(ns,condsD,'^-', 'LineWidth', 2,'DisplayName',"Chafee-Infante", "MarkerSize",10)

load("data/data_icesheet","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'*-', 'LineWidth', 2,'DisplayName',"ice sheet", "MarkerSize",10)
figure(2)
semilogy(ns,condsD,'*-', 'LineWidth', 2,'DisplayName',"ice sheet", "MarkerSize",10)

load("data/data_burgers","O_errors","condsD");
ns = 1:numel(O_errors);
figure(1)
semilogy(ns,O_errors,'+-', 'LineWidth', 2,'DisplayName',"Burgers'", "MarkerSize",10)
figure(2)
semilogy(ns,condsD,'+-', 'LineWidth', 2,'DisplayName',"Burgers'", "MarkerSize",10)

figure(1)
ylabel("relative operator error","Interpreter","latex", "FontSize",15)
xlabel("ROM dimension","Interpreter","latex","FontSize",15)
set(gca, 'YScale', 'log')
grid on
legend("show","Interpreter","latex","FontSize",12)

savefig("figures/operator_errors.fig")
exportgraphics(gcf,"figures/operator_errors.pdf")

figure(2)
ylabel("condition number","Interpreter","latex","FontSize",15)
xlabel("ROM dimension","Interpreter","latex","FontSize",15)
set(gca, 'YScale', 'log')
grid on
legend("show","Interpreter","latex","FontSize",12)

savefig("figures/condition_numbers.fig")
exportgraphics(gcf,"figures/condition_numbers.pdf")
