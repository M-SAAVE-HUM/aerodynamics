%% M-SAAVE Preliminary Sizing Plot
% n = W/S = wing loading
n_vector = linspace(0,100,1001);

% Environmental Parameters
rho = 1.225;    % air density (kg/m^3)
g = 9.81;       % gravity (m/s^2)

% Aircraft Parameters
mass = 25;      % mass (kg)
W = mass * g;   % weight (N)
Sref = 1.5;     % wing area (m^2)
AR = 6;         % aspect ratio 
e = 0.95;       % efficiency factor
CL_max = 1.4;   % max CL of aircraft (from AVL)
CD0 = 0.05;     % parasitic drag
T = 130;        % max static thrust (N)

%% Requirements
% Stall Requirement
v_stall = sqrt((2*W) / (rho*CL_max*Sref)); % stall speed (m/s)
n_stall = 0.5 * rho * (v_stall^2) * CL_max;

% Takeoff Requirement
k_to = 1.2; % ratio of takeoff speed to stall speed
takeoff_dist = 100; % takeoff distance requirement (self set, m)
takeoff_TW = n_vector .* (k_to^2 / g*rho*takeoff_dist*CL_max);

% Landing Requirement

% Climb Requirement
k_c = 1.2; % ratio of climb speed to stall speed
K = 
G = 
climb_TW = (k_c^2 / CL_max)*CD0 + (CL_max / k_c^2)*K + G;

% Cruise Requirement

% Ceiling Requirement

% Maneuver Requirement

%% T/W - W/S Plot (Thrust-to-Weight vs. Wing Loading)
hold on
plot(n_vector, takeoff_TW)
%plot(n_vector, turn_TW)
%plot(n_vector, cruise_TW)
yline(climb_TW, 'm')
xline(n_stall, 'g')
scatter(W/Sref, T/W)
%xlim([0 100])
%ylim([0 .5])
%legend("Takeoff", "Turn", "Cruise", "Climb", "Stall", "Design");
xlabel('Wing Loading (W/S)') % figure out units
ylabel('Thrust-to-Weight Ratio (T/W)') % figure out units
hold off