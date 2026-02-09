function [ac] = aircraft_1()
% Sample aircraft configuration file

%% Misc

% Relative locations of components (to wing root LE)
ac.x = 0.8;  % Distance between wing and hstab AC's [m]
ac.xv = 0.8; % Distance between wing and vstab AC's [m]


%% Mass Properties

ac.mass = 10.7;  % Total mass of aircraft [kg]

% CG location
ac.xcg = 0.1419; % CG x-coordinate relative to LE of root [m]
ac.ycg = 0.0; % CG y-coordinate relative to LE of root [m]
ac.zcg = 0.0; % CG z-coordinate relative to LE of root [m]

% Inertia matrix [kg/m^3]
% ac.Ixx = ;
% ac.Iyy = ;
% ac.Izz = ;
% ac.Ixy = ;
% ac.Iyz = ;
% ac.Ixz = ;


%% Wing

% Planform
wing.Sref = 0.8013;     % Reference planform area [m^2]
wing.AR = 7.87;         % Aspect ratio [-]
wing.taper1 = 1.0;      % Inboard taper ratio [-]
wing.taper2 = 0.5;      % Outboard taper ratio [-]
wing.btaper = 0.5;      % b/2 fraction of kink [-]
wing.dihedral1 = 0.0;   % Inboard dihedral angle [deg]
wing.dihedral2 = 10.0;  % Outboard dihedral angle [deg]
wing.incidence = 0.0;   % Wing incidence angle [deg]

wing.bail_st = 0.5;     % Aileron section start b/2 fraction [-]
wing.bail_end = 0.9;    % Aileron section end b/2 fraction [-]
wing.ca_c = 0.22;       % Aileron chord fraction [-]

wing.xLE = 0.0;         % Wing root leading edge x-location [m]
wing.yLE = 0.0;         % Wing root leading edge y-location [m]
wing.zLE = 0.0;         % Wing root leading edge z-location [m]

% Airfoils
% Root
wing.airfoil_root = load('model_airplane\airfoils\e201.mat');   % Load airfoil data from .mat file
wing.aroot = wing.airfoil_root.filename;        % Name of airfoil .dat file
wing.CLAF_root = wing.airfoil_root.CLAF;    % Cl_alpha modifier factor
wing.CDCL_root = wing.airfoil_root.CDCL;    % Sectional polar points
% Kink
wing.airfoil_kink = wing.airfoil_root;
wing.akink = wing.airfoil_kink.filename;
wing.CLAF_kink = wing.airfoil_kink.CLAF;
wing.CDCL_kink = wing.airfoil_kink.CDCL;
% Aileron start
wing.airfoil_ail_st = wing.airfoil_root;
wing.aail_st = wing.airfoil_ail_st.filename;
wing.CLAF_ail_st = wing.airfoil_ail_st.CLAF;
wing.CDCL_ail_st = wing.airfoil_ail_st.CDCL;
% Aileron end
wing.airfoil_ail_end = wing.airfoil_root;
wing.aail_end = wing.airfoil_ail_end.filename;
wing.CLAF_ail_end = wing.airfoil_ail_end.CLAF;
wing.CDCL_ail_end = wing.airfoil_ail_end.CDCL;
% Tip
wing.airfoil_tip = wing.airfoil_root;
wing.atip = wing.airfoil_tip.filename;
wing.CLAF_tip = wing.airfoil_tip.CLAF;
wing.CDCL_tip = wing.airfoil_tip.CDCL;

% DEPENDENTS - Calculate
wing.b = sqrt(wing.AR*wing.Sref);      % Wingspan [m]
wing.b1 = wing.btaper*wing.b/2;        % Root-kink length [m]
wing.b2 = wing.b/2 - wing.b1;          % Kink-tip length [m]
wing.croot = wing.Sref/(wing.b1*(1+wing.taper1) + wing.b2*wing.taper1*(1+wing.taper2));  % Root chord length [m]
wing.ckink = wing.taper1*wing.croot;            % Chord length at kink [m]
wing.ctip = wing.taper2*wing.taper1*wing.croot; % Tip chord length [m]

cbar1 = 2/3*wing.croot*(1 + wing.taper1 + wing.taper1^2)/(1+wing.taper1);  % MAC length for first section [m]
S1 = wing.b1*(wing.croot+wing.ckink);                                      % Reference area of first section [m^2]
cbar2 = 2/3*wing.ckink*(1 + wing.taper2 + wing.taper2^2)/(1+wing.taper2);  % MAC length for second section [m]
S2 = wing.b2*(wing.ckink+wing.ctip);                                       % Reference area of second section [m^2]
wing.cbar = (cbar1*S1 + cbar2*S2)/wing.Sref;                               % MAC length for whole wing [m]
wing.taperav = (wing.taper1*wing.btaper + wing.taper2*(1 - wing.btaper))/(wing.taper1 + wing.taper2);   % Weighted average taper ratio [-]
wing.sweep_LE = atand(tand(0) - 4/wing.AR*((0 - 0.25)/100*(1-wing.taperav)/(1+wing.taperav)));          % Average leading edge sweep [deg]
wing.x_cbar = (wing.b/6)*(wing.croot+2*wing.ctip)/(wing.croot+wing.ctip)*tand(wing.sweep_LE);           % x-location of MAC leading edge relative to root LE [m]


% Store in ac
ac.wing = wing;


%% Horizontal-stabilizer

% Planform
hstab.AR = 6;           % Aspect ratio [-]
hstab.taper = 0.5;      % Taper ratio [-]
hstab.dihedral = 0.0;   % Dihedral angle [-]
hstab.incidence = 0.0;  % H-stab incidence angle [-]
hstab.yLE = 0.0;        % y-location of root leading edge [m]
hstab.zLE = 0.05;        % z-loaction of root leading edge [m]

% Elevator
hstab.ce_crh = 0.4;     % Elevator chord fraction [-] at root

% Airfoil
hstab.airfoil = load('model_airplane\airfoils\NACA0012.mat'); % Load airfoil data   
hstab.CLAF = hstab.airfoil.CLAF;      % Cl_alpha modifier factor
hstab.CDCL = hstab.airfoil.CDCL;      % Sectional polar points


% DEPENDENTS - Calculate
cht = 0.70; % General aviation - single engine (CAN CHANGE) [-]
Lht = ac.x; % Lever arm from wing AC to hstab AC [m]
hstab.Sref = cht*wing.cbar*wing.Sref/Lht;   % H-stab Sref [m^2]

hstab.b = sqrt(hstab.AR*hstab.Sref);      % Wingspan of hstab [m]
hstab.croot = 2*hstab.Sref/(hstab.b*(1 + hstab.taper)); % Root chord [m]
hstab.ctip = hstab.croot*hstab.taper;     % Tip chord [m]

hstab.sweep_LE = atand(tand(0) - 4/hstab.AR*((0 - 0.25)/100*(1-hstab.taper)/(1+hstab.taper)));       % Leading edge sweep angle [deg]
hstab.cbar = (2/3)*(hstab.croot + hstab.ctip - (hstab.croot*hstab.ctip)/(hstab.croot + hstab.ctip)); % Mean Aerodynamic Chord [m]
hstab.x_cbar = (hstab.b/6)*(hstab.croot+2*hstab.ctip)/(hstab.croot+hstab.ctip)*tand(hstab.sweep_LE); % x-location of LE of MAC relative to root LE [m]
hstab.xLE = ac.x + wing.xLE + wing.x_cbar - hstab.x_cbar;   % x-location of root LE [m]


% Store in ac
ac.hstab = hstab;



%% Vstab

% Planform
vstab.AR = 4;       % Aspect ratio [m^2]
vstab.taper = 0.6;  % Taper ratio [-]
vstab.yLE = 0.0;    % y-coord of LE [m]
vstab.zLE = 0.0;    % z-coord of LE [m]

% Rudder
vstab.ce_crv = 0.4; % Rudder chord fraction at root [-]

% Airfoil
vstab.airfoil = load('model_airplane\airfoils\NACA0012.mat'); % Load airfoil data
vstab.CLAF = vstab.airfoil.CLAF;      % Cl_alpha modifier factor
vstab.CDCL = vstab.airfoil.CDCL;      % Sectional polar points


% DEPENDENTS - Calculate
cvt = 0.04; % General aviation - single engine (CAN CHANGE) [-]
Lvt = ac.xv; % Lever arm from wing AC to hstab AC
vstab.Sref = cvt*wing.b*wing.Sref/Lvt;

vstab.b = sqrt(vstab.AR*vstab.Sref);
vstab.croot = 2*vstab.Sref/(vstab.b*(1 + vstab.taper));
vstab.ctip = vstab.croot*vstab.taper; 

vstab.sweep_LE = atand(tand(0) - 4/vstab.AR*((0 - 0.25)/100*(1-vstab.taper)/(1+vstab.taper)));
vstab.cbar = (2/3)*(vstab.croot + vstab.ctip - (vstab.croot*vstab.ctip)/(vstab.croot + vstab.ctip));
vstab.x_cbar = (vstab.b/6)*(vstab.croot+2*vstab.ctip)/(vstab.croot+vstab.ctip)*tand(vstab.sweep_LE);
vstab.xLE = ac.xv + wing.xLE + wing.x_cbar - vstab.x_cbar;

ac.vstab = vstab;

end