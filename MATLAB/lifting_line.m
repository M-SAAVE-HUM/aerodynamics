%% Lifting Line Calculator
% uses lifting line theory to calculate drag polar for given geometry
clear; close all;

% inputs
% AR
% airfoil
% airfoil stall angle
% speed
% reynolds number
airfoil_type = 'naca';          % set airfoil type (naca or coords)
number = num2str(6412);         % set the airfoil number
%number = load('ch10.txt');     % OR import coordinates
alpha = deg2rad(12);            % cruise angle of attack
Uinf = 22;                      % cruise speed in m/s
AR = 6;                         % aspect ratio

% outputs
% L vs y, cl vs y, cd vs cl plots
% CL and L, CD and D at 0 aoa
% CL max and angle
% Max L/D and angle

% along the way
% lifting line at cruise and many other angles
% estimates for parasitic drag

disp(repmat('-', 1, 60));
disp('Results!')
disp(repmat('-', 1, 60));
disp('At aoa = °, CL = , L = N')
disp('At aoa = °, CD = , D = N')
disp('CLmax = , occurs at aoa = °')
disp('Max L/D = , occurs at aoa = °')
disp(repmat('-', 1, 60));