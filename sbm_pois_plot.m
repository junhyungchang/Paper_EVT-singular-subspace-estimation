load("sbm_pois.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
histogram(targetsig,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
title("Poisson SBM")
legend("$T_{\mathbf{s},n}$", "Gumbel PDF", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
text_x = 5.5; 
text_y = 0.35;
text(text_x, text_y, '$p=\log n/\sqrt{n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
% text(text_x, text_y, '$p=\sqrt{\log n/n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
% text(text_x, text_y, '$p=1/\sqrt{n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
text_x = 5.25;
text_y = 0.3;
text(text_x, text_y, '$q=p/\log n$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
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
xlabel('Standard Gumbel Quantiles', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

histogram(targetc,40, 'Normalization','pdf', 'FaceColor','k', 'FaceAlpha', 0.3)
hold on
x = linspace(-4,8,1e+3);
v = 1;
pdf = v*exp(-v*x-exp(-v*x));
plot(x,pdf, 'Color', 'k', 'LineWidth', 1.2)
title("Poisson SBM")
legend("$\widetilde{T}_{\tilde{\mathbf{s}},n}$", "Gumbel pdf", "Interpreter","latex", "FontSize", 20)
% title("$n=2e+3$, $m=5e+2$, $r=10$",'Interpreter','latex')
xlim([-4,8])
ylim([0,0.5])
text_x = 5.5; 
text_y = 0.34;
text(text_x, text_y, '$p=\log n/\sqrt{n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
% text(text_x, text_y, '$p=\sqrt{\log n/n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
% text(text_x, text_y, '$p=1/\sqrt{n}$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
text_x = 5.25;
text_y = 0.29;
text(text_x, text_y, '$q=p/\log n$', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter','latex');
set(gca,"FontSize",14)
hold off

figure(6)
sortdatac = sort(targetc);
plot(theoretical_quantiles, sortdatac, 'k+', 'MarkerFaceColor', 'b');
hold on;
min_val = min([theoretical_quantiles; sortdatac]);
max_val = max([theoretical_quantiles; sortdatac]);
plot([min_val, max_val], [min_val, max_val], '--', 'Color', [0.5, 0.5, 0.5],  'LineWidth', 2);
xlabel('Standard Gumbel Quantiles', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);
xlim([-3,8])
ylim([-3,8])
set(gca,"FontSize",14)
hold off