%% Aircraft Performance!!!
% thrust, power, range, endurance
% from kfid 201 notes

% aircraft parameters
g = 9.81; % gravity (m/s^2)
rho = 1.225; % air density (kg/m^3)
Re = 10^6; % Reynolds number
mass = 25; % kg 
W = mass*g;
S = 1.5; % wing area (m^2)
AR = 6.5; % aspect ratio
e = 0.85; % oswald efficiency factor
CL_max = 1.4; % maximum CL (@ stall, alpha = 12 degrees) 

% velocities (m/s)
v_stall = sqrt(W / (0.5 * rho * S * CL_max)); % stall speed (m/s)
v_cruise = 24; % cruise speed (m/s)
v = linspace(v_stall,40,101); 

% drag estimation
[CD0, CDi_est, CD_est, D_est] = drag_estimator(rho,v_cruise,Re,e,AR,S,mass); 
CD_stall = 0.15817;

% thrust wrt velocity
T = (0.5.*rho.*v.^2.*S*CD0) + (2.*W^2)./(rho.*v.^2.*S*pi*AR*e);
test_T = 134 - 3.68.*v + 0.0275.*v.^2;

% using dynamic thrust test values, find optimal cruise speed
syms u
drag_curve = @(u) .5.*CD_stall.*rho.*(u^2).*S;
T_curve = @(u) 134 - 3.68.*u + 0.0275.*u.^2; % dynamic thrust curve
v_opt_cruise = fzero(@(u) T_curve(u) - drag_curve(u), 1);

% power wrt velocity - check the efficiency (n) values
n_prop = 0.8; 
n_motor = 0.9;
eta = n_prop*n_motor;
P = (T.*v)./eta;

% min thrust
V_min_thrust = ((4*W^2) / (rho^2 * S^2 * pi * AR * e * CD0))^0.25;
min_thrust = min(T);

% minimum power 
P_min = (4/3) * sqrt((2*W^3) / (rho*S)) * ((3*CD0) / (pi*AR*e)^3)^0.25;
[min_power, p_index] = min(P);
V_min_power = v(p_index);

% plots
tcl = tiledlayout(1,2);
s = strcat('Aircraft Performance');
title(tcl, s)

nexttile
hold on;
plot(v, T, '-r', 'Linewidth', 1.5);
xline(v_stall,'-k','Stall Speed');
xline(v_cruise,'-k','Cruise Speed');
title('Thrust vs. Velocity')
ylabel('Thrust (N)');
xlabel('Velocity (m/s)');
grid on;
hold off;

nexttile
plot(v, P, '-b', 'Linewidth', 1.5);
xline(v_stall,'-k','Stall Speed');
xline(v_cruise,'-k','Cruise Speed');
title('Power vs. Velocity')
ylabel('Power (W)');
xlabel('Velocity (m/s)');
grid on;

X = sprintf('Minimum thrust = %s N, at v = %s m/s', ...
    num2str(min_thrust), num2str(V_min_thrust));
disp(X);

Y = sprintf('Minimum power = %s W, at v = %s m/s', ...
    num2str(min_power), num2str(V_min_power));
disp(Y);
