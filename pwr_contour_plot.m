% Create power contour plot (Figure 3)
r = 10;
n_len = 20;
nmax =400; nmin =20;
load("pwr_contour.mat") 
n = powermat(:,1);

% two-to-infinity norm effect size
d_n = powermat(:,2);
d_n(d_n < 1e-15) = 0;

% Z values (empirical power)
z = powermat(:,3);

% Define the grid resolution
% num_n_points = 11;
num_d_n_points = 40;

% Create the grid
n_lin = linspace(20, 400, 20);
d_n_lin = linspace(0, max(d_n), num_d_n_points);
[Nq, Dq] = meshgrid(n_lin, d_n_lin);

% Interpolate the data onto the grid
Zq = griddata(n, d_n, z, Nq, Dq, 'linear');

locind = Dq > 2* sqrt(r./Nq);
Zq(locind == 1) = NaN;

% Calculate contour levels and round them to two decimal places
num_levels = 4; % Number of contour levels
Z_min = min(z);
Z_max = max(z);
levels = linspace(Z_min, Z_max, num_levels);
levels = round(levels * 100) / 100; % Round to two decimal places
% Remove duplicate levels (if any) after rounding
levels = unique(levels);

% Compute the theoretical threshold curve d_n = 2 * sqrt(r / n)
n_plot = linspace(nmin,nmax,100);
d_n_plot = 2 * sqrt(r ./ n_plot);
% lam1 = 2./(n_plot.^(1/2).*log(n_plot).^(2/3)).^2; % weak signal
lam1 = 2./(n_plot).^2; % strong signal
% lam1 = 2./(n_plot.^(3/2));
a_n = sqrt(lam1)./(2*sqrt(log(2*n_plot)));
% b_n = sqrt(lam1.*log(2*n_plot))-(log(log(2*n_plot))+log(4*pi)).*a_n/2;
b_n = sqrt(lam1.*log(2*n_plot));


%fill in the empty areas
cmap = colormap;
numLevels = size(cmap, 1);
colorIdx = round((1-min(Zq(:))) / (max(Zq(:)) - min(Zq(:))) * (numLevels - 1)) -10 ; % 245
fillColor = cmap(colorIdx, :);
% fillColor = cmap(100, :);
fill([n_plot, fliplr(n_plot)], [d_n_plot, zeros(size(d_n_plot))], fillColor, 'EdgeColor', 'none');
hold on;

% Create the filled contour plot without contour lines
figure(1);
[~, h_filled] = contourf(Nq, Dq, Zq, 20, 'LineColor', 'none');

% Adjust colormap to lighter colors to enhance contour line visibility
cmap = colormap('jet');
cmap_lighter = 0.6 + 0.4 * cmap; % Scale colors to be lighter
colormap(cmap_lighter);

% Overlay contour lines indicating Z values
[contour_matrix, h_lines] = contour(Nq, Dq, Zq, levels, 'LineColor', 'k', 'LineWidth', 0.5);

% Label the contour lines with Z values
clabel(contour_matrix, h_lines, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 200);

% actually plot the theoretical lines
plotub = plot(n_plot, d_n_plot, 'k:', 'LineWidth', 1.7);
% plotan = plot(n_plot, a_n/sqrt(1.4*r), 'k-.', 'LineWidth', 2);
plotbn = plot(n_plot, b_n*sqrt(1.4*r)*0.25, 'k--', 'LineWidth', 2);
legend([plotub, plotbn], {'$2\sqrt{r\mu/n}$', '$0.25b_n (\mu r)^{1/2}$'}, 'Interpreter', 'latex', 'Location', 'northeast');

% Set labels and title
xlabel('$n$', 'Interpreter','latex');
ylabel('$d_n$', 'Interpreter','latex');
title('Rejection rates under $\mathrm{H}_\mathrm{A}$ (strong signal)', "Interpreter","latex"); 
% % <- strong signal
% title('Rejection rates under $\mathrm{H}_\mathrm{A}$ (weak signal)', "Interpreter","latex");
% % <- weak signal

% Adjust plot limits if necessary
xlim([min(n), max(n)]);
ylim([0, 0.5]);

% % Optional: Add colorbar
colorbar;

% Ensure that the background of the plot is white
set(gca, 'Color', 'white', 'FontSize', 14);
hold off;
