% General script to creat an AVL model (geometry and mass) for a given 
% configuration

% REQUIRED:
%   - An airfoil polar file exported from XFOIL (as .txt) that contains 
%     negative stall, positive stall, 0 deg alpha, and 5 deg alpha. This
%     file should be in the 'airfoils' folder.
%   - An aircraft configuration file (.m) with all configuration data
%     put in. This should be in the 'airplanes' folder.

% ASSUMPTIONS:
%   - 0 quarter-chord sweep (unneccesary an difficult to manufacture for a 
%     subsonic UAV)
%   - 0 geometric twist (also hard to effectively manufacture)
%   - Constant chord fraction ailerons (so stringers and/or spars are
%     straight, also common practice)
%   - Straight elevator and rudder hingelines (not constant chord
%     fraction) (easier manufacturing and allows for one servo to be used).
%     This can cause issues if the taper and/or root chord fraction are
%     incompatable.
%   - Can contain one 'kink' or 'crank' in the planform i.e. can be double
%     tapered/dihedraled
%   - The aileron must be at the kink or outboard of it
%   - No flaps (unneccesary?)
%   - No fuselage (yet)
%   - Point mass with no rotational inertia (yet)
%   - The horizontal stabilizer must not be inline with the wing planform
%     i.e. it must have a different z-translation (trailing horseshoe
%     vortice legs can run into the h-stab and cause computational errors)
%   - The parasitic drag estimate from the 'CDCL' AVL function is a rough
%     ballpark to produce trends and the actual value itself should be
%     viewed with heavy skepticism.
%   - XFOIL can be optimistic about Clmax/min's and under-predict Cd, if
%     wind tunnel data is available for relevant Re's, they should be used 
%     instead.
%   - No propeller interaction modeled (duh, AVL can't do that)
%   - If 'NACA' is in the airfoil name, AVL will use its own geometry
%     insetad of a .dat file.

% RECOMMENDTIONS FOR IMPROVEMENT
%   - Make wing planform input a vector of parameters for the sections so
%     more kinks/cranks can be produced.
%   - Make a function to execute AVL and produce a drag polar and detect
%     stall using method of sections.
%   - Make a robust airfoil analysis function and place it inside the 
%     aircraft configuration file to automatically run XFOIL and produce 
%     reliable data.
%     - Perhaps even run at different Re's depending on the chord length.
%   - Use this within fmincon() to roughly optimize a configuration.
%   - Better folder/file management so loading into AVL is cleaner.

% CREATED BY:
%     Andrew Painter, 2/5/26

% UPDATED BY:
%     Andrew Painter, 2/9/26

clear
clc
close all

name = "test";

%% Run airfoils

file = 'model_airplane/airfoils/e201.txt'; % Pacc file
af = read_pacc(file); % Save into .mat file

file = 'model_airplane/airfoils/NACA0012.txt'; % Pacc file
af = read_pacc(file); % Save into .mat file


%% Produce aircraft

ac = aircraft_1(); % Generate full aircraft configuration struct


%% Produce geometry file

geom_folder = "model_airplane/geometries"; % Destination folder
gen_plane(ac,name,geom_folder); % Generate .avl file


%% Produce mass file

rho = 1.225;    % Air density (can be changed later) [kg/m^3]
mass_folder = "model_airplane/masses"; % Destination folder
gen_mass(ac,rho,name,mass_folder); % Generate .mass file

fprintf("\n%s configuration modeled, have fun with loading it!\n",name);


%% AVL load example

% load model_airplane/geometries/test.avl
% mass model_airplane/masses/test.mass
