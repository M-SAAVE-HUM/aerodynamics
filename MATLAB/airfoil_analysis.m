%% Airfoil Analysis

% airfoil parameters
m = mfoil('naca','2412', 'npanel',199);

%% optional mfoil plots
% basic plot
% m.setoper('alpha', 10, 'Re', 8e5, 'Ma', 0.07);
% m.solve

% set a related easier operating condition  
m.setoper('alpha', 4, 'Re', 8e5);

% run the solver  
m.solve

% turn off BL initialization 
m.oper.initbl = false;

% set the difficult operating condition  
% m.setoper('alpha',14.5, 'Re',10^6);
% 
% % run the solver: now converges (right, bottom) 
% m.solve 

% plot surface variable distributions at chosen aoa
%m.plot_distributions;

% plot streamlines
%m.plot_stream([-.1, 1.1,-.2,.2],51);

%% analysis over more conditions
% alphas = linspace(-5, 10, 101);
% CLs = zeros(1, length(alphas));
% CDs = zeros(1, length(alphas));
% CMs = zeros(1, length(alphas));
% for i = 1:length(alphas)
%     m.setoper('alpha', alphas(i), 'Re', 8e5, 'Ma', 0.07);
%     m.solve
%     CLs(i) = m.post.cl;
%     CDs(i) = m.post.cd;
%     CMs(i) = m.post.cm;
% end
% 
% % plot drag polar
% figure;
% plot(CDs, CLs, 'LineWidth', 2);
% xlabel('$C_d$', 'fontsize', 14, 'interpreter', 'latex')
% ylabel('$C_l$', 'fontsize', 14, 'interpreter', 'latex')
% 
% % plot cl vs alpha
% figure;
% plot(alphas, CLs, 'LineWidth', 2);
% xlabel('angle of attack $\alpha$', 'fontsize', 14, 'interpreter', 'latex')
% ylabel('$C_l$', 'fontsize', 14, 'interpreter', 'latex')
% 
% % plot cm vs alpha
% figure;
% plot(alphas, CMs, 'LineWidth', 2);
% xlabel('angle of attack $\alpha$', 'fontsize', 14, 'interpreter', 'latex')
% ylabel('$C_m$', 'fontsize', 14, 'interpreter', 'latex')
% 
% % output cl, cd, cm at specific aoas
% % 0
% % cruise
% % climb