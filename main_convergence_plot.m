load("main_convergence.mat")
% Change name if .mat file name is changed in main_convergence_mat.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the three histogram-qqplot pairs in Figure 1.
figure(1)
histogram(targetsig,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$T_{\mathbf{s},n}$ (Oracle)", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(4)
pd = makedist('GeneralizedExtremeValue',0,1,0);
sortdata = sort(targetsig);
m = length(sortdata);
ps = (1:m)' / (m + 1); % Empirical cumulative probabilities
theoretical_quantiles = icdf(pd, ps);
plot(theoretical_quantiles, sortdata, 'k+');
hold on;
min_val = min([theoretical_quantiles; sortdata]);
max_val = max([theoretical_quantiles; sortdata]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
histogram(targetuc,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$\widetilde{T}_{\hat{\mathbf{s}},n}^{\mathrm{uc}}$ (Uncorrected)", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(5)
sortdatauc = sort(targetuc);
plot(theoretical_quantiles, sortdatauc, 'k+', 'MarkerFaceColor', 'b');
hold on;
min_val = min([theoretical_quantiles; sortdatauc]);
max_val = max([theoretical_quantiles; sortdatauc]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

histogram(real(targetc),40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
legend("$\widetilde{T}_{\tilde{\mathbf{s}},n}$ (De-biased)", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
set(gca,"FontSize",14)
hold off

figure(6)
sortdatac = sort(real(targetc));
plot(theoretical_quantiles, sortdatac, 'k+', 'MarkerFaceColor', 'b');
hold on;
min_val = min([theoretical_quantiles; sortdatac]);
max_val = max([theoretical_quantiles; sortdatac]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Theoretical Quantiles (Standard Gumbel)', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off