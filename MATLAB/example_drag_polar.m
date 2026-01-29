%% Example Drag Polar
% Kabir Khwaja, 2/27/25

clear; clc; close all;

% Aerodynamic data (made up values)
CD_data = [0.023 0.021 0.024 0.029 0.035 0.045 0.05 0.06 0.07 0.08 0.09];
CL_data = [0.3 0.45 0.6 0.8 1.0 1.15 1.3 1.45 1.55 1.6 1.65];

% Fit a quadratic function CD = a*CL^2 + b*CL + c
p = polyfit(CL_data, CD_data, 2);

% Generate model data for plotting
CL_model = linspace(-0.5, 1.75, 100); % can change these bounds if needed
CD_model = polyval(p, CL_model);

% Compute L/D ratio and find max L/D
LD_ratio = CL_model ./ CD_model;
[LD_max, idx_max] = max(LD_ratio);
CD_LDmax = CD_model(idx_max);
CL_LDmax = CL_model(idx_max);

% Compute tangent line at max L/D
slope = LD_max;
CD_tangent = linspace(0, CD_LDmax*1.5, 100);
CL_tangent = slope * CD_tangent;

% Estimate CD_min and corresponding CL
[CD_min, idx_CDmin] = min(CD_model);
CL_minD = CL_model(idx_CDmin);

% Make plot!
figure;
hold on;

% Plot experimental data
scatter(CD_data, CL_data, 75,'xk', 'LineWidth', 1.25, 'DisplayName', 'Data');

% Plot drag model
plot(CD_model, CL_model, '-k', 'LineWidth', 1.5, 'DisplayName', 'Drag Model');

% Plot (L/D)_max point
scatter(CD_LDmax, CL_LDmax, 75, 'or', 'filled', 'DisplayName', 'L/D_{max}');

% Plot tangent line
plot(CD_tangent, CL_tangent, '--r', 'LineWidth', 1.5, 'DisplayName', 'Tangent');

% Highlight (CD_min, CL_minD)
scatter(CD_min, CL_minD, 75, 'ok', 'filled');

% Labels and legend
xlabel('$C_D$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$C_L$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Data', 'Drag Polar', 'Max L/D', 'Tangent', '$C_{Dmin}$', 'Location', 'northwest', 'Interpreter', 'latex');

hold off;

% Save figure as EPS file
saveas(gcf, 'drag_polar_plot.eps', 'epsc');
