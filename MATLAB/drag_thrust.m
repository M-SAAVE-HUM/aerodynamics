%% Plot drag vs velocity
close all; clear; clc;
m = 25; % mass, kg
AR = 6; % aspect ratio
S = 1.5; % wing area, m^2
e = 0.9; % efficiency factor
rho = 1.225; % air density, kg/m^3
mu = 1.81*10^-5; % viscosity, Pa*s
L = 0.5; % use chord length as characteristic length

v = linspace(10, 40, 61);
CD0 = zeros(1, length(v));
CDi = zeros(1, length(v));
CD = zeros(1, length(v));
D = zeros(1, length(v));

for i = 1:length(v)
    v_curr = v(i);
    Re = rho*v_curr*L/mu;
    [CD0_curr, CDi_curr, CD_curr, D_curr] = drag_estimator(rho, v_curr, Re, e, AR, S, m);
    CD0(i) = CD0_curr;
    CDi(i) = CDi_curr;
    CD(i) = CD_curr;
    D(i) = D_curr;
end

%% Dynamic thrust curves
% thrust data
Tspeeds = 0.44704*[10, 20, 30, 40, 50, 60]; % wind tunnel speeds converted to m/s
T20 = 2.*[2.3421776, -1.2114684, -3.081021468, -4.394366807, -6.271936219, -8.57836442];
T30 = 2.*[10.83466887, 6.820875415, 1.408542487, -1.591695877, -5.48043084, -9.20726663];
T40 = 2.*[20.89059396, 17.2615482, 11.74422025, 6.812580539, 1.62399816, -3.930862795];
T50 = 2.*[35.38376927, 31.29135204, 26.47852589, 21.43823048, 13.65938878, 6.0330272];
T60 = 2.*[49.39413752, 44.60040068, 39.81684575, 34.2094657, 26.68805023, 18.20731627];
T70 = 2.*[62.51138192, 56.98849228, 52.19456227, 45.85557424, 38.6756476, 30.02118989];
T80 = 2.*[73.9846495, 68.94531098, 64.85434904, 57.32135527, 49.76234755, 41.39913283];

% thrust curve polyfitting
T20_curve = polyfit(Tspeeds,T20,2);
T30_curve = polyfit(Tspeeds,T30,2);
T40_curve = polyfit(Tspeeds,T40,2);
T50_curve = polyfit(Tspeeds,T50,2);
T60_curve = polyfit(Tspeeds,T60,2);
T70_curve = polyfit(Tspeeds,T70,2);
T80_curve = polyfit(Tspeeds,T80,2);

% create thrust curve functions for plotting
f20 = polyval(T20_curve,v);
f30 = polyval(T30_curve,v);
f40 = polyval(T40_curve,v);
f50 = polyval(T50_curve,v);
f60 = polyval(T60_curve,v);
f70 = polyval(T70_curve,v);
f80 = polyval(T80_curve,v);

%% TODO: drag curves for cruise, climb (get e factor from AVL)

figure(1)
hold on;
plot(v, CD0, 'linewidth', 1.5);
plot(v, CDi, 'linewidth', 1.5);
plot(v, CD, 'linewidth', 1.5);
xlabel('Velocity [m/s]', 'fontsize', 14, 'interpreter', 'latex')
ylabel('Drag Quantities', 'fontsize', 14, 'interpreter', 'latex')
legend('Parasitic Drag Coefficient, $C_{D0}$', 'Induced Drag Coefficient, $C_{Di}$', 'Drag Coefficient, $C_D$', 'fontsize', 12, 'interpreter', 'latex', 'location', 'best')
hold off;

figure(2)
hold on;
plot(v, D, 'linewidth', 1.5);
plot(v, CD0.*(0.5*rho.*v.^2), 'linewidth', 1.5);
plot(v, CDi.*(0.5*rho.*v.^2), 'linewidth', 1.5);
xlabel('Velocity [m/s]', 'fontsize', 14, 'interpreter', 'latex')
ylabel('Forces [N]', 'fontsize', 14, 'interpreter', 'latex')
legend('Total Drag', 'Parasitic Drag', 'Induced Drag', 'fontsize', 12, 'interpreter', 'latex', 'location', 'best')

% drag curve with all thrust curves
figure(3)
hold on;
plot(v, D, 'linewidth', 1.5);
plot(v, f20, 'linewidth', 1.5);
plot(v, f30, 'linewidth', 1.5);
plot(v, f40, 'linewidth', 1.5);
plot(v, f50, 'linewidth', 1.5);
plot(v, f60, 'linewidth', 1.5);
plot(v, f70, 'linewidth', 1.5);
plot(v, f80, 'linewidth', 1.5);
xlabel('Velocity [m/s]', 'fontsize', 14, 'interpreter', 'latex')
ylabel('Forces [N]', 'fontsize', 14, 'interpreter', 'latex')
legend('Total Drag', 'Thrust (20\% throttle)', 'Thrust (30\% throttle)', ...
    'Thrust (40\% throttle)', 'Thrust (50\% throttle)', ...
    'Thrust (60\% throttle)', 'Thrust (70\% throttle)', ...
    'Thrust (80\% throttle)', ...
    'fontsize', 12, 'interpreter', 'latex', 'location', 'best')

% drag curve with just 60, 70, 80% throttle curves
figure(4)
hold on;
plot(v, D, 'linewidth', 1.5);
plot(v, f50, 'linewidth', 1.5);
plot(v, f60, 'linewidth', 1.5);
plot(v, f70, 'linewidth', 1.5);
plot(v, f80, 'linewidth', 1.5);
xlabel('Velocity [m/s]', 'fontsize', 14, 'interpreter', 'latex')
ylabel('Forces [N]', 'fontsize', 14, 'interpreter', 'latex')
legend('Total Drag', 'Thrust (50\% throttle)', 'Thrust (60\% throttle)', 'Thrust (70\% throttle)', ...
    'Thrust (80\% throttle)', 'fontsize', 12, 'interpreter', 'latex', 'location', 'best')