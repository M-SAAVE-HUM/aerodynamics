%% drag polar
clear; close all;

alphas = -5:1:12;
CLs = [0.0719, 0.1533, 0.2346, 0.3159, 0.397, 0.4779, 0.5585, 0.6388, ...
    0.7187, 0.7981, 0.877, 0.9553, 1.033, 1.11, 1.1862, 1.2617, 1.3362, 1.41];
CDs = [0.03383, 0.03439, 0.03586, 0.03768, 0.04041, 0.04385, 0.048, 0.05285, ...
    0.05836, 0.06459, 0.07146, 0.07897, 0.0871, 0.09584, 0.10514, 0.115, 0.1254, 0.13627];

plot(CDs, CLs, 'LineWidth', 2)
xlabel('$C_D$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$C_L$', 'FontSize', 20, 'Interpreter', 'latex')
%title('Drag Polar (Payload Configuration)', 'FontSize', 20, 'Interpreter', 'latex')

% Compute CL/CD
E = CLs ./ CDs;

% Find maximum efficiency
[Emax, idx] = max(E);

% Get corresponding values
CD_maxE = CDs(idx);
CL_maxE = CLs(idx);
alpha_maxE = alphas(idx);

% Plot the point
hold on;
plot(CD_maxE, CL_maxE, 'ro', 'MarkerSize', 10, 'LineWidth', 2)

% Annotate the point
text(CD_maxE+0.003, CL_maxE-0.01, sprintf('Max L/D', ...
    Emax, alpha_maxE), ...
    'FontSize', 14, 'Interpreter', 'latex', 'Color', 'r');
hold off;

figure(2)
L_D = CLs./CDs;
plot(alphas, L_D, 'LineWidth', 2)
xlabel('$\alpha$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$L/D$', 'FontSize', 20, 'Interpreter', 'latex')
disp(alphas(10))
hold on;
plot(4, 12.36, 'ro', 'MarkerSize', 10, 'LineWidth', 2)

% Annotate the point
text(4+0.4, 12+0.7, sprintf('Max L/D', ...
    Emax, alpha_maxE), ...
    'FontSize', 14, 'Interpreter', 'latex', 'Color', 'r');
hold off;