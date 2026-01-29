%% Takeoff Distance Calculator!!!
% made with love by performance team <3

% from Raymer:
% takeoff dist = ground roll + rotation up + transition + climb
% TD = SG + SR + STR + SC

clear; close all

% environment parameters
g = 9.81;       % gravity [m/s^2]
rho = 1.225;    % air density [kg/m^3]
Re = 800000;    % Reynolds number

% aircraft parameters
altitude = 100;                 % altitude (m)
mass = 25;                      % mass of aircraft at takeoff (kg)
mass_landing = 20;              % mass at landing (no payload, kg)
W = mass*g;                     % weight at takeoff (N)
W_landing = mass_landing*g;     % weight at landing (N)
S = 1.5;                        % wing area (m^2)
b = 3;                          % span (m)
h = 0.25;                       % wing height above ground (m)
AR = b^2 / S;                   % aspect ratio
e = 0.85;                       % efficiency factor (CHECK THIS VALUE)
CL_max = 1.5;                   % max CL of craft (CHECK THIS VALUE)

% velocity bounds
v_stall = sqrt(W / (0.5 * rho * S * CL_max)); % stall speed (m/s)
v_takeoff = 1.1 * v_stall; % takeoff speed
v_landing = 1.15 * v_stall; % landing speed
v_cruise = 23; % m/s
CL = (2*W) / (rho * S * v_takeoff^2);

% [coeff parasitic drag, coeff induced drag, coeff drag, total drag force]
[CD0, CDi, CD, D] = drag_estimator(rho, v_stall, Re, e, AR, S, mass);

%% Method 1: Integration
syms v % make sure you have symbolic toolbox installed

% ground lift and drag coefficients
CLg = CL * 0.9; 
CDg = CD * 1.1;

Lg = 0.5*rho*S*CLg*v^2;  % lift force
Dg = 0.5*rho*S*CDg*v^2;  % drag force

% Dynamic thrust curve (80% throttle)
aT = -0.0545; 
bT = -1.2093;
cT = 154.1938;
T = aT*v^2 + bT*v + cT;

% Friction
mu = 0.08;
Ff = mu*(W - Lg);

% acceleration
accel = (T - Dg - Ff)/mass;
Sg_sym = int(v/accel, v, 0, v_takeoff);
SG_int = double(vpa(Sg_sym, 3));

% Time to takeoff
T_takeoff_sym = int(1/accel, v, 0, v_takeoff);
T_takeoff = double(vpa(T_takeoff_sym, 3));

%% Method 2: Raymer
% simpler formula, more of an estimate
e = (1 + b^2/(256*h^2))*e; % modified for ground effect 
K = 1/(pi*AR*e); % induced drag factor
T0 = 120; % takeoff thrust [N]
KT = T0/W - mu; 
KA = (rho*S)/(2*W) * (mu*CL_max - CD0 - K*CL_max^2);
SG = (1/(2*g*KA)) * log((KT + KA*v_takeoff^2)/KT); % ground roll

%% Add up all phases of takeoff
SR = v_takeoff; % rotate up, raymer: assume takes 1 second

v_transition = 1.15 * v_stall;
n = 1.2; % load factor
R = (v_transition^2) / (0.2*g); % radius of transition
sin_gamma = (T0-D)/W; % sine of angle created during transition
gamma = asin(sin_gamma);
hTR = R*(1-cos(gamma)); % height of transition
STR = R*sin_gamma; % transition distance

h_obstacle = 100; % obstacle height (m)
SC = (h_obstacle - hTR)/tan(gamma);

% if aircraft clears obstacle during transition use this, take out climb
%STR = sqrt(R^2 - (R-h_obstacle)^2); % transition distance

% times for diff phases
total_takeoff = SG + SR + STR + SC;
T_SG = T_takeoff; % ground roll time estimate
T_SR = 1; % raymer assumption = 1 sec
T_STR = STR/v_transition;
T_SC = SC/(v_cruise);
total_time = T_SG + T_SR + T_STR + T_SC;

%% abracadabra
X = sprintf('Method 1 (integration): ground roll distance = %s m', ...
    num2str(SG_int));
disp(X);
disp(['Time to reach takeoff speed = ', num2str(T_takeoff), ' s']);

Y = sprintf('Method 2 (raymer): ground roll distance = %s m', ...
    num2str(SG));
disp(Y);

Z = sprintf('total takeoff distance = %s m', ...
    num2str(total_takeoff));
disp(Z);
disp(['Total takeoff time (under assumptions) = ', num2str(total_time), ' s']);

extra = sprintf(['takeoff distance by segments: ' ...
    'SG = %s m, SR = %s m, STR = %s m, SC = %s m'], ...
    num2str(SG), num2str(SR), num2str(STR), num2str(SC));
disp(extra);

extra2 = sprintf(['takeoff time by segments: ' ...
    'SG = %s s, SR = %s s, STR = %s s, SC = %s s'], ...
    num2str(T_SG), num2str(T_SR), num2str(T_STR), num2str(T_SC));
disp(extra2);
