function [] = gen_wing(wing,name,template,folder)
% DESCRIPTION:
%   Generates an AVL geometry file for a double-tapered wing.
%   Assumes: no sweep, no geometric twist.
%
% INPUTS:
%   wing -> struct of wing parameters (see below) [struct]
%   name -> name of .avl file (e.g. "wing_1") [string]
%   template -> relative file location and name of geometry template file
%   folder -> relative folder location for output file
%
% OUTPUT:
%   [VOID] -> produces .avl file in specified folder
%
% CREATED:
%   Andrew Painter, 2/3/26
%
% LAST MODIFIED:
%   2/3/26


%% Unload variables

template = fileread(template);  % Load template

Sref = wing.Sref;       % Reference planform area [m^2]
AR = wing.AR;           % Aspect Ratio [-]
taper1 = wing.taper1;   % Taper ratio 1 (c_root/c_mid-station) [-]
taper2 = wing.taper2;   % Taper ratio 2 (c_mid-station/c_tip) [-]
btaper = wing.btaper;   % Half-span fraction of mid-station
dihedral1  = wing.dihedral1;    % Dihedral of root-mid-span section [deg]
dihedral2  = wing.dihedral2;    % Dihedral of mid-span-tip section [deg]
incidence = wing.incidence;     % Angle of incidence of wing relative to fuselage centerline [deg]
xLE = wing.xLE;         % Root leading edge x-location
yLE = wing.yLE;         % Root leading edge y-location
zLE = wing.zLE;         % Root leading edge z-location

% Root airfoil and aero model
aroot = wing.aroot;         % Root airfoil filename
CLAF_root = wing.CLAF_root; % Root airfoil lift-slope factor
CDCL_root = wing.CDCL_root;

% Kink airfoil and aero model
akink = wing.akink;         % Kink airfoil filename
CLAF_kink = wing.CLAF_kink; % Root airfoil lift-slope factor
CDCL_kink = wing.CDCL_kink;


% Tip airfoil and aero model
atip  = wing.atip;          % Tip airfoil filename
CLAF_tip = wing.CLAF_tip;   % Root airfoil drag polar
CDCL_tip = wing.CDCL_tip;

    

%% Further definitions

b = sqrt(AR*Sref);      % Wingspan [m]
b1 = btaper*b/2;        % Root-mid-span length [m]
b2 = b/2 - b1;          % Mid-span-tip length [m]
croot = Sref/(b1*(1+taper1) + b2*taper1*(1+taper2));  % Root chord length [m]
ckink = taper1*croot;            % Chord length at mid-span kink [m]
ctip = taper2*taper1*croot;     % Tip chord length [m]


cbar1 = 2/3*croot*(1 + taper1 + taper1^2)/(1+taper1);   % MAC length for first section [m]
S1 = b1*btaper*croot*(1+taper1);                        % Reference area of first section [m^2]
cbar2 = 2/3*ckink*(1 + taper2 + taper2^2)/(1+taper2);    % MAC length for second section [m]
S2 = b2*croot*(1+taper1);                               % Reference area of second section [m^2]
cbar = (cbar1*S1 + cbar2*S2)/Sref;                      % MAC length for whole wing [m]

Xref = xLE + croot*0.25;        % Reference x-coordinate for moment calculation [m]

dXw  = xLE;     % x-translate wing surface
dYw  = yLE;     % y-translate wing surface
dZw  = zLE;     % z-translate wing surface

%% Define coordinates

% Root
yroot = yLE;
xroot = xLE;
zroot = zLE;
troot = 0.0;

% Kink
ykink = b1;
xkink = xLE + 0.25*(croot - ckink);
zkink = zLE + b1*sind(dihedral1);
tkink  = 0.0;

% Tip
ytip = b/2;
xtip = xLE + 0.25*(croot - ctip);
ztip = zkink + b2*sind(dihedral2);
ttip = 0.0;


%% Map replacements

replacements = {

    "{{NAME}}",          name;
    "{{Sref}}",          sprintf("%.6f", Sref);
    "{{cbar}}",          sprintf("%.6f", cbar);
    "{{bref}}",          sprintf("%.6f", b);
    "{{Xref}}",          sprintf("%.6f", Xref);

    "{{dXw}}",           sprintf("%.6f", dXw);
    "{{dYw}}",           sprintf("%.6f", dYw);
    "{{dZw}}",           sprintf("%.6f", dZw);
    "{{incidence}}",     sprintf("%.6f", incidence);

    "{{xroot}}",         sprintf("%.6f", xroot);
    "{{yroot}}",         sprintf("%.6f", yroot);
    "{{zroot}}",         sprintf("%.6f", zroot);
    "{{croot}}",         sprintf("%.6f", croot);
    "{{troot}}",         sprintf("%.6f", troot);
    "{{aroot}}",         sprintf(aroot);
    "{{CDCL_root}}",     sprintf(CDCL_root);
    "{{CLAF_root}}",     sprintf("%.6f", CLAF_root);

    "{{xkink}}",         sprintf("%.6f", xkink);
    "{{ykink}}",         sprintf("%.6f", ykink);
    "{{zkink}}",         sprintf("%.6f", zkink);
    "{{ckink}}",         sprintf("%.6f", ckink);
    "{{tkink}}",         sprintf("%.6f", tkink);
    "{{akink}}",         sprintf(akink);
    "{{CDCL_kink}}",     sprintf(CDCL_kink);
    "{{CLAF_kink}}",     sprintf("%.6f", CLAF_kink);

    "{{xtip}}",          sprintf("%.6f", xtip);
    "{{ytip}}",          sprintf("%.6f", ytip);
    "{{ztip}}",          sprintf("%.6f", ztip);
    "{{ctip}}",          sprintf("%.6f", ctip);
    "{{ttip}}",          sprintf("%.6f", ttip);
    "{{atip}}",          sprintf(atip);
    "{{CDCL_tip}}",      sprintf(CDCL_tip);
    "{{CLAF_tip}}",      sprintf("%.6f", CLAF_tip);

};


%% Replace template

output = template;
for k = 1:size(replacements,1)
    output = strrep(output, replacements{k,1}, replacements{k,2});
end


%% Output file

filename = name + ".avl";      % Output filename
directory = filename;
fileID = fopen(directory,'w');
fwrite(fileID,output);
fclose(fileID);

fprintf("AVL geometry file generated to: %s\n",filename);


end