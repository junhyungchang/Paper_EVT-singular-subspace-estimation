load("t_dist.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create three histogram + QQ plot pairs for Figure 4
figure(1)
histogram(target1,300, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$t_4$ noise", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(4)
pd = makedist('GeneralizedExtremeValue',0,1,0);
sortdata = sort(target1);
m = length(sortdata);
ps = (1:m)' / (m + 1); % Empirical cumulative probabilities
theoretical_quantiles = icdf(pd, ps);
plot(theoretical_quantiles, sortdata, 'k+');
hold on;
min_val = min([theoretical_quantiles; sortdata]);
max_val = max([theoretical_quantiles; 7.5]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,11])
set(gca,"FontSize",14)
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
histogram(target2,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$t_5$ noise", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(5)
sortdatauc = sort(target2);
plot(theoretical_quantiles, sortdatauc, 'k+', 'MarkerFaceColor', 'b');
hold on;
min_val = min([theoretical_quantiles; sortdatauc]);
max_val = max([theoretical_quantiles; 7.5]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

histogram(target3 ,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$t_{10}$ noise", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(6)
sortdatac = sort(target3);
plot(theoretical_quantiles, sortdatac, 'k+', 'MarkerFaceColor', 'b');
hold on;
min_val = min([theoretical_quantiles; sortdatac]);
max_val = max([theoretical_quantiles; 7.5]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off