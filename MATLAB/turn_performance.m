%% Turning Performance!!!
% made with love by Kabir, 3/6/24

% given craft values and turning angle, find turn radius and time
% two methods: estimate with kinematics, estimate using bank angle
% both should be approx equal for bank angle = 45 degrees
% bank angle estimate generally more accurate!!!

clear
close all

% parameters
rho = 1.225; % air density (kg/m^3)
CLmax = 1.25; % full aircraft CLmax
m = 25; % mass of aircraft (kg)
g = 9.81;
S = 1.5; % wing reference area (m^2)
u = 18; % turning speed (m/s)
t_ang = 90; % turn angle (degrees)
b_ang = 30; % bank angle (degrees)

[r1, t1, r2, t2, tr1, tr2] = turning(rho, m, u, S, CLmax, t_ang, b_ang);

% abracadabra
A = sprintf('turning radius = %s m, turning time = %s sec, rate = %s rad/s' ...
    , num2str(r1), num2str(t1), num2str(tr1));
disp(A);

% abracadabra 2 electric boogaloo
B = sprintf('alt: turning radius = %s m, turning time = %s sec, rate = %s rad/s' ...
    , num2str(r2), num2str(t2), num2str(tr2));
disp(B);

%% Turning Function
% input angles in degrees
function[tr1, tt1, tr2, tt2, w, psi] = turning(rho, m, u, S, CL, ta, ba)

    % setup
    g = 9.81; % gravity
    W = m*g; % weight of craft
    t_angle = deg2rad(ta); % turning angle
    phi = deg2rad(ba); % bank angle
    
    % Method 1
    % calculations based on kinematics
    L_max_turn = 0.5 * CL * rho * u^2 * S; % max lift in turn (N)
    Fc = sqrt( L_max_turn^2 - W^2 ); % centripetal force (N)
    tr1 = (m * u^2)/Fc; % turning radius (m)
    ac = u^2/tr1; % centripetal acceleration (m/s^2)
    w = sqrt(ac/tr1); % turning rate (rad/s)
    tt1 = t_angle/w; % turning time (s)

    % Method 2
    % calculations based on bank angle - adjusts turning speed
    L = W/cos(phi);
    n = L/W; % load factor
    v = sqrt( (n*W) / (0.5 * rho * S * CL) ); % turning speed 
    tr2 = (1 / (0.5 * rho * g * CL) ) * (W/S) * (n / sqrt(n^2 - 1));
    psi = (g/v) * sqrt(n^2 - 1);
    tt2 = t_angle/psi;

end