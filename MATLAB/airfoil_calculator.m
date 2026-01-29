%% Airfoil Calculator!!!    
% made with love by Kabir <3 (Winter '24)

% This code is for preliminary airfoil analysis. It provides a plot of the
% airfoil selected, and performs basic analysis of a rectangular wing with
% the selected airfoil. You can input the aspect ratio and span,
% and the program will output plots of the airfoil, drag polar, lift
% distribution across the wing, lift coefficient vs angle of attack, and 
% lift to drag ratio vs angle of attack. The code can also be modified to
% find the span required to produce a given lift force.

% This program uses mfoil to determine values relevant to chosen airfoil, 
% so can also be used for airfoils outside NACA series from coordinates.
% Documentation at: https://public.websites.umich.edu/~kfid/codes.html 

%% sections:
% 1: Input definitions and parameters
% 2: Lifting Line calculation at cruise
% 3: Crunch numbers from parts 2 and 3 to get cruise lift and drag
% 4: Wetted area approximation for parasitic drag
% 5: Repeat parts 2 and 3 for different angles of attack
% 6: Process data, make pretty plots :)

clear
close all

%% 1: Definitions and Parameters
airfoil_type = 'naca';          % set airfoil type (naca or coords)
number = num2str(6412);         % set the airfoil number
%number = load('ch10.txt');     % OR import coordinates
alpha = deg2rad(12);            % cruise angle of attack
Uinf = 22;                      % cruise speed in m/s
AR = 8;                         % aspect ratio

% input a desired Lift to determine span required
g = 9.81;                       % gravity (m/s^2)
mass = 25;                      % mass of aircraft (kg)
L_des = mass*g;                 % desired lift at cruise = weight

% OR, input known span to determine amount of lift produced
b = 3;                          % span = 3m
S = b^2 / AR;                   % wing area (m^2)

% VERY IMPORTANT!!! find the zero lift angle of attack
% can do this using a separate mfoil script (zero_lift.m)
alphaL0 = deg2rad(-6.2907);

rho = 1.225; % air density (kg/m^3)
gamma = 1.4; % specific heat ratio
R = 287.185; % specific gas constant (J/kgK)
T = 288.15; % temperature (K)
Re = 750000; % Reynold's number
Mach = Uinf / sqrt(gamma*R*T); % Mach number

%% 2: kfid <3
m = mfoil(airfoil_type, number, 'npanel',199); % load mfoil
m.setoper('Re', Re)
m.setoper('Ma', Mach)
% run the solver
m.solve 

% Prandtl's Lifting Line Theory yurrr
a0 = 2*pi; % 2D lift curve slope - assumption for naca airfoils
Nu = 20; % number of unknowns
M = zeros(Nu,Nu); % coefficient matrix
F = zeros(Nu,1); % right−hand−side vector

% span locations at which to enforce Prandtl’s equation
theta = linspace(pi/2, pi-pi/(2*Nu), Nu); % set up system
for i = 1:Nu % loop over spanwise locations
% loop over unknowns
    for j = 1:Nu
        M(i,j) = 4*AR/a0*sin((2*j-1)*theta(i)) ...
        + (2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i)); 
    end
    % fill in right−hand side
    F(i) = alpha - alphaL0;
end

% Solve the system - will use results in section 4
A = M\F;

%% 3: crunching numbers 
% Lift coefficient
A1 = A(1); % fourier series coefficients!!!
A2 = A(2);
A3 = A(3);
CL = pi*AR*A1;
L = 0.5*rho*S*CL*Uinf^2; % lift

% Induced drag coefficient
delta = 0; 
for j = 2:Nu
    delta = delta + (2*j-1)*(A(j)/A1)^2; 
end
e = 1/(1+delta); % span efficiency factor
CDi = CL^2/(pi*AR)*(1+delta); % induced drag

% find span for wing to produce desired lift at cruise speed
% b = sqrt((2*L) / (pi * rho * A1 * Uinf^2)) % span (meters) yippeeee!!!
% c = span/AR % chord (meters)

% point collocation sectional lift distribution
theta = linspace(0, pi, 100);
Y = -b/2 * cos(theta);
lift = (2*b*rho*Uinf^2) .* (A1.*sin(theta) + A3.*sin(3*theta));
cl_dist = lift./(0.5*rho*0.5*Uinf^2);
CL_dist = (((pi*AR)./sin(theta)) .* (A1.*sin(theta) + A3.*sin(3*theta)));

%% 4: aerodynamics is such a drag :(
% total drag = parasitic + induced
% parasitic drag calculated using wetted area approximation, uses drag
% estimator function

%[CD0, CDi_est, CD_est, D_est] = drag_estimator(rho,Uinf,Re,e,AR,S,mass); 
CD0 = 0.0286;
CD = CDi + CD0; % total aircraft drag

%% 5: i love data!!!
% get CL and CD values for many angles of attack, use for drag polar
K = 31; % number of iterations
alphas = linspace(-15, 15, K); % plot between -15 and 15 degrees AoA
co_lift = zeros(size(alphas));
co_drag = zeros(size(alphas));
lift_drag_ratio = zeros(size(alphas));
for k = 1:K
    m = mfoil(airfoil_type, number, 'npanel',199); % load mfoil
    % set angle of attack
    m.setoper('alpha',deg2rad(alphas(k)));
    m.setoper('Ma', Mach)
    % disabled Reynolds number so iterating doesn't take forever
    %m.setoper('Re',Re)
    % run the solver, plot the results 
    m.solve, clf; 

    % repeat PLL
    M_new = zeros(Nu,Nu); 
    F_new = zeros(Nu,1); 
    
    % span locations
    theta = linspace(pi/2, pi-pi/(2*Nu), Nu); 
    for i = 1:Nu 
    % loop over unknowns
        for j = 1:Nu
            M_new(i,j) = 4*AR/a0*sin((2*j-1)*theta(i)) ...
            + (2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i)); 
        end
        % fill in right−hand side
        F_new(i) = deg2rad(alphas(k)) - alphaL0; 
    end

    % Solve the new system
    A_new = M_new\F_new;
    co_lift(k) = pi*AR*A_new(1);
    delta = 0;
    for j = 2:Nu
        delta = delta + (2*j-1)*(A_new(j)/A_new(1))^2; 
    end
    co_drag(k) = (co_lift(k).^2) ./ (pi*AR)*(1+delta) + CD0;
    lift_drag_ratio(k) = co_lift(k)/co_drag(k);

end    

% calculate max CL
[max_lift, max_lift_index] = max(co_lift);
max_lift_alpha = alphas(max_lift_index);

% calculate CLminD
% lift coefficient at minimum drag
[min_value, min_index] = min(co_drag);
CLminD = co_lift(min_index);

% calculate max L/D
[max_LD, max_LD_index] = max(lift_drag_ratio);
max_LD_alpha = alphas(max_LD_index);

%% 6: post processing
% another full CD craft estimate
K = 1 / (pi*AR*e);
CD_full = CD0 + K*(CL - CLminD)^2; % total drag coefficient
D = 0.5*rho*S*CD_full*Uinf^2; % total drag force

% mfoil plot of airfoil at input conditions
m = mfoil(airfoil_type, number, 'npanel',199); 
% if desired, replace alpha below with max_alpha
m.setoper('alpha', alpha); 
m.setoper('Ma', Mach)
m.setoper('Re', Re)
nexttile
m.solve 

% sectional lift vs y - plot vs. elliptical
fig2 = figure(2);
y_dist = linspace(-1.5, 1.5, 3001);
cl_chord_max_ellipse = 2 * 0.57 * 1.5 / (pi * 1.5);
elliptical = cl_chord_max_ellipse * sqrt(1 - (y_dist ./ 1.5).^2) .* (0.5*rho*Uinf^2);

hold on;
plot(Y, lift, 'b-', 'Linewidth', 2);
plot(y_dist, elliptical, 'k--', 'Linewidth', 2);
xlabel('Span location, $y$ [$m$]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Sectional Lift [$N/m$]', 'Interpreter', 'latex', 'FontSize', 16);
legend('Initial Design Wing', 'Elliptical', 'Interpreter', 'latex', 'Location', 'northeast')
hold off;
saveas(fig2, 'Lprime.eps', 'epsc');

% sectional lift coefficient vs y
fig3 = figure(3);
hold on;
plot(Y, cl_dist, 'r-', 'Linewidth', 2);
yline(1.65, 'k--', 'Linewidth', 2)
xlabel('Span location, $y$ [$m$]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Sectional Lift Coefficient, $C_l$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Spanwise $C_l$', 'NACA 6412 $C_{l,max}$', 'Interpreter', 'latex', 'Location', 'northeast')
hold off;
saveas(fig3, 'cl_dist.eps', 'epsc');

% drag polar
fig4 = figure(4);
plot(co_drag, co_lift, 'b', 'Linewidth', 2);
xlabel('Drag Coefficient, $C_D$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Lift Coefficient, $C_L$', 'Interpreter', 'latex', 'FontSize', 16);
saveas(fig4, 'drag_polar.eps', 'epsc');

% reports
W = sprintf('At alpha = %s, CL = %s, Lift = %s N', ...
    num2str(alpha), num2str(CL), num2str(L));
disp(W);

X = sprintf('At alpha = %s, CD estimate 1: %s, CD estimate 2: %s', ...
    num2str(alpha), num2str(CD), num2str(CD_full));
disp(X);

Y = sprintf('CLmax = %s, occurs at alpha = %s degrees', ...
    num2str(max_lift), num2str(max_lift_alpha));
disp(Y);

Z = sprintf('Max L/D = %s, occurs at %s degrees', ...
    num2str(max_LD), num2str(max_LD_alpha));
disp(Z)
