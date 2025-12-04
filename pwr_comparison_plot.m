% Parameters
r = 5; % Example value for r (adjust as needed)
n_len = 20;
nmax =400; nmin =20;
load("pwr_comparison.mat","powerourmat","powerfrobmat");
powerourmat = sortrows(powerourmat, [-1,1]); % change the row orders for inverted plot
powerfrobmat = sortrows(powerfrobmat, [-1,1]); 

% row discrepency parameter
mu = 400-powerfrobmat(:,1);
mu1 = 400-powerourmat(:,1);

% signal strength vector
s = powerfrobmat(:,2);
s1 = powerourmat(:,2);
ss = unique(s);

% Z values (empirical power)
z = powerfrobmat(:,3);
z1 = powerourmat(:,3);
z2 = powerourmat(:,3) - powerfrobmat(:,3);

muunique = unique(mu);
sunique = unique(s);
[X,Y] = meshgrid(sunique,muunique);
Z = reshape(z2,length(sunique), length(muunique))';
muunique1 = unique(mu1);
sunique1 = unique(s1);
[X1,Y1] = meshgrid(sunique1,muunique1);
Z1 = reshape(z1,length(sunique1), length(muunique1))';

% Calculate contour levels and round them to two decimal places
num_levels = 5; % Number of contour levels
% Need to change Z_min and Z_max for each plot
% Z_min = min(z);  % for individual plot (legacy)
% Z_max = max(z);
Z_min = min(z2);  % for difference plot
Z_max = max(z2);
levels = linspace(Z_min, Z_max, num_levels);
levels = round(levels * 100) / 100; % Round to two decimal places
levels = [0, levels]; % for difference plot
% levels = [0.2, 0.6, 1]; % for individual plot
% Remove duplicate levels (if any) after rounding
levels = unique(levels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the filled contour plot without contour lines
figure(1);
[~, h_filled] = contourf(X, Y, Z, 20, 'LineColor', 'none');

% Adjust colormap to lighter colors to enhance contour line visibility
cmap = colormap('jet');
cmap_lighter = 0.6 + 0.4 * cmap; % Scale colors to be lighter
colormap(cmap_lighter);

% % Adjust transparency for regions with no data
% set(h_filled, 'AlphaData', ~isnan(Zq));
% set(h_filled, 'FaceAlpha', 'flat');
% set(h_filled, 'AlphaDataMapping', 'none');

% Overlay contour lines indicating Z values
hold on;
[contour_matrix, h_lines] = contour(X, Y, Z, levels, 'LineColor', 'k', 'LineWidth', 0.5);

% Label the contour lines with Z values
clabel(contour_matrix, h_lines, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 200);


% Set labels and title
xlabel('$s_{r}$', 'Interpreter','latex');
ylabel('$\xi$', 'Interpreter','latex');
title('Difference in rejection rates of $\widetilde{T}_{\tilde{\mathbf{s}},n}$ and $\widetilde{T}_{\mathrm{Frob}}$', "Interpreter","latex");

% Adjust plot limits
xlim([log10(400^(0.5)), max(s)]);
ylim([min(mu), max(mu)]);

% Define the vector of x-tick values
xTickValues = (ss(3:2:end)); 
yTickValues = linspace(0,20,6);
y_minor_ticks = 0:1:20;
% yTickValues(1) = 1;
% Define the corresponding labels
xTickLabels = ["$n^{0.55}$", "$n^{0.7}$", "$n^{0.85}$", "$n$"];

% Set ticks
xticks(xTickValues);
yticks(yTickValues);
% Set x-axis tick labels
xticklabels(xTickLabels);

% Enable LaTeX interpreter for labels
ax = gca;
ax.YMinorTick = 'on';
ax.YAxis.MinorTickValues = y_minor_ticks;
ax.TickLength = [0.02, 0.01];
ax.TickLabelInterpreter = 'latex';

% add vertical line at SNR threshold
n = 400;
% vline = xline(log10(sqrt(n)*log(n)), '--k', 'LineWidth', 2);
% legend(vline, "$s = \sqrt{n\log(n)}$", "Interpreter","latex", "Location","southeast")

% % Optional: Add colorbar
colorbar;

% Ensure that the background of the plot is white
set(gca, 'Color', 'white', 'FontSize', 15);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the filled contour plot without contour lines
figure(2);
[~, h_filled] = contourf(X1, Y1, Z1, 20, 'LineColor', 'none');

% Adjust colormap to lighter colors to enhance contour line visibility
cmap = colormap('jet');
cmap_lighter = 0.6 + 0.4 * cmap; % Scale colors to be lighter
colormap(cmap_lighter);

% % Adjust transparency for regions with no data
% set(h_filled, 'AlphaData', ~isnan(Zq));
% set(h_filled, 'FaceAlpha', 'flat');
% set(h_filled, 'AlphaDataMapping', 'none');

% Overlay contour lines indicating Z values
hold on;
[contour_matrix, h_lines] = contour(X1, Y1, Z1, levels, 'LineColor', 'k', 'LineWidth', 0.5);

% Label the contour lines with Z values
clabel(contour_matrix, h_lines, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 200);


% Set labels and title
xlabel('$s_{r}$', 'Interpreter','latex');
ylabel('$\xi$', 'Interpreter','latex');
title('Rejection rates for $\widetilde{T}_{\tilde{\mathbf{s}},n}$', "Interpreter","latex");

% Adjust plot limits if necessary
xlim([log10(400^(0.5)), max(s1)]);
ylim([min(mu1), max(mu1)]);

% Define the vector of x-tick values
xTickValues = (ss(3:2:end)); 
yTickValues = linspace(0,20,6);
% yTickValues(1) = 1;
% Define the corresponding labels
xTickLabels = ["$n^{0.55}$", "$n^{0.7}$", "$n^{0.85}$", "$n$"];

% Set ticks
xticks(xTickValues);
yticks(yTickValues);
% Set x-axis tick labels
xticklabels(xTickLabels);

% Enable LaTeX interpreter for labels
ax = gca;
ax.YMinorTick = 'on';
ax.YAxis.MinorTickValues = y_minor_ticks;
ax.TickLength = [0.02, 0.01];
ax.TickLabelInterpreter = 'latex';

% add vertical line at SNR threshold
n = 400;
% vline = xline(log10(sqrt(n)*log(n)), '--k', 'LineWidth', 2);
% legend(vline, "$s = \sqrt{n\log(n)}$", "Interpreter","latex", "Location","southeast")

% % Optional: Add colorbar
colorbar;

% Ensure that the background of the plot is white
set(gca, 'Color', 'white', 'FontSize', 15);
hold off;
