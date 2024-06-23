%% Takeoff Distance Calculator!!!
% made with love by performance team <3

% from Raymer:
% takeoff dist = ground roll + rotation up + transition + climb
% TD = SG + SR + STR + SC

clear
close all

% environment parameters
g = 9.81; % GRAVITY (m/s^2)
rho = 1.225; % air density (kg/m^3)
Re = 10^6; 

% aircraft parameters
altitude = 525; % altitude (m)
mass = 25; % mass of aircraft at takeoff (kg)
mass_landing = 20; % mass at landing (no payload, kg)
W = mass*g; % weight at takeoff (N)
W_landing = mass_landing*g; % weight at landing (N)
S = 1.5; % wing area (m^2)
b = 3; % span (m)
h = 0.25; % wing height above ground (m)
AR = b^2 / S; % aspect ratio
CL_max = 1.4; % max CL of craft

% velocity bounds
v_stall = sqrt(W / (0.5 * rho * S * CL_max)); % stall speed (m/s)
v_takeoff = 1.1 * v_stall; % takeoff speed
v_landing = 1.15 * v_stall; % landing speed
CL = (2*W) / (rho * S * v_takeoff^2);

% [coeff parasitic drag, coeff induced drag, coeff drag, total drag force]
[CD0, CDi, CD, D] = drag_estimator(rho, v_stall, 10^6, 0.85,6.5,1.5,mass);

%% Method 1: Integration
syms v % make sure you have symbolic toolbox installed

% ground lift and drag coefficients
CLg = CL * 0.9; % lift slightly lower
CDg = CD * 1.1; % drag slightly higher
Lg = 0.5*rho*S*CLg*v^2; % lift
Dg = 0.5*rho*S*CDg*v^2; % drag

% thrust data (change this to fit dynamic thrust data)
T0 = 130; % static thrust 100% (N)
T = T0 * 0.75; % thrust at takeoff (estimate)

% Dynamic thrust data
% a = -0.0344; 
% b = -1.0;
% c = 50.9;
% T = a*v^2 + b*v + c; % Total Thrust

e_init = 0.85; % efficiency factor
e = (1 + b^2/(256*h^2))*e_init; % modified for ground effect
K = 1/(pi*AR*e); % induced drag factor

% Friction
mu = 0.05; % friction coefficient - raymer
Ff = mu*(W-Lg); % Total friction force

% integrate
accel = (g/W)*(T-Dg-Ff);
Sg = int(v/accel, v, 0, v_takeoff);
Sg = vpa(Sg, 3);
SG_int = double(Sg);

%% Method 2: Raymer
% simpler formula, more of an estimate
KT = T/W - mu; 
KA = (rho*S)/(2*W) * (mu*CL - CD0 - K*CL^2);
SG = (1/(2*g*KA)) * log((KT + KA*v_takeoff^2)/KT); % ground roll

%% Add up all phases of takeoff
SR = v_takeoff; % rotate up, raymer: assume takes 1 second

v_transition = 1.15 * v_stall;
n = 1.2; % load factor
R = (v_transition^2) / (0.2*g); % radius of transition
sin_gamma = (T-D)/W; % sine of angle created during transition
gamma = asin(sin_gamma);
hTR = R*(1-cos(gamma)); % height of transition
STR = R*sin_gamma; % transition distance

h_obstacle = 20; % obstacle height (m)
SC = (h_obstacle - hTR)/tan(gamma);

% if aircraft clears obstacle during transition use this, take out climb
%STR = sqrt(R^2 - (R-h_obstacle)^2); % transition distance

total_takeoff = SG + SR + STR + SC;

%% abracadabra
X = sprintf('method 1 (integration): ground roll distance = %s m', ...
    num2str(SG_int));
disp(X);

Y = sprintf('method 2 (raymer): ground roll distance = %s m', ...
    num2str(SG));
disp(Y);

Z = sprintf('total takeoff distance = %s m', ...
    num2str(total_takeoff));
disp(Z);

extra = sprintf(['takeoff distance by segments: ' ...
    'SG = %s m, SR = %s m, STR = %s m, SC = %s m'], ...
    num2str(SG), num2str(SR), num2str(STR), num2str(SC));
disp(extra);
