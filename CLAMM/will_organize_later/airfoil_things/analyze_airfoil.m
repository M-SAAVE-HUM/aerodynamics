function [CDCL,CLAF,data] = analyze_airfoil(airfoil,Re,M,plot)
% DESCRIPTION:
%   Produces a polar for an airfoil using mfoil to be used in the aero
%   model in AVL.
%   Requires mfoil.m
%
% INPUTS:
%   airfoil -> Airfoil .dat file string OR NACA code (e.g. 'NACA 2412') [string]
%   Re -> Reynolds number [-]
%   M  -> Mach number [-]
%   plot (optional) -> Plot drag polar toggle [1 -> plotting on, 0 -> plotting off]
%
% OUTPUTS:
%   CDCL -> Vector of drag polar values for use in AVL [1x3]
%   CLAF -> Lift-curve slope factor for use in AVL [1x1]
%   data -> Detailed result struct [struct]
%
% CREATED:
%   Andrew Painter, 2/3/26
%
% MODIFIED:
%   Andrew Painter, 2/3/26
%

arguments
    airfoil string
    Re double
    M  double
    plot (1,1) {mustBeNumeric} = 0; 
end


%% Run mfoil
clear m

% Setup result coefficient vectors
Cls = [];   % Sectional lift coefficient [-]
Cds = [];   % Sectional drag coefficient [-]
Cms = [];   % Sectional moment coefficient [-]


% Load airfoil
airfoil_des = char(airfoil);
if contains(airfoil_des,'NACA') == 1
    NACA = strtrim(erase(airfoil,'NACA'));
    try 
        m = mfoil('naca',NACA, 'npanel',199);
    catch E
        disp(E)
        fprintf("Error: %s\n",E.message);
        return
    end

else
    af = load(airfoil);
    airfoil = strtrim(erase(airfoil,'.dat'));
    try 
        m = mfoil('coords', af, 'npanel',199);
    catch E
        printf("Error: %s\n",E.message);
        return
    end
end

% Set settings
m.param.verb = 0;      % Set output verbosity (0=quiet, 3=verbose)
m.param.niglob = 100;   % Set number of Newton iterations (50 is default)
m.param.rtol = 1e-6;   % Set convergence tolerance (1e-10 is default)

alphas = [];

% Initialize BL solution with 0.1 deg alpha
m.oper.initbl = true;
m.setoper('alpha',0.25, 'Re',Re, 'Ma',M);
m.solve;
m.oper.initbl = false;

% Run positive side of polar
alphas_pos = [0.5:0.25:10 10.1:0.1:20];
for i=1:length(alphas_pos)
    alphas_pos(i)
    m.setoper('alpha',alphas_pos(i), 'Re',Re, 'Ma',M);
    try
        lastwarn('');
        m.solve;
        [msg, ~] = lastwarn;
        
        if ~isempty(msg)
            fprintf("Warning detected at %.2f deg: %s\n", alpha, msg);
            break
        end
    catch Esolve
        fprintf("Failure detected, sweep stopped.\n")
        break;
    end
    Cls = [Cls, m.post.cl];
    Cds = [Cds, m.post.cd];
    Cms = [Cms, m.post.cm];
    alphas = [alphas, alphas_pos(i)];
end


% Initialize BL solution with 0 deg alpha
m.oper.initbl = true;
m.setoper('alpha',0, 'Re',Re, 'Ma',M);
m.solve;
m.oper.initbl = false;
Cls = [Cls, m.post.cl];
Cds = [Cds, m.post.cd];
Cms = [Cms, m.post.cm];

% Run negative side of alphas
alphas_neg = [-0.25:-0.25:-10 -10.1:-0.1:-20];
for i=1:length(alphas_neg)
    alphas_neg(i)
    m.setoper('alpha',alphas_neg(i), 'Re',Re, 'Ma',M);
    try
        lastwarn('');
        m.solve;
        [msg, ~] = lastwarn;
        
        if ~isempty(msg)
            fprintf("Warning detected at %.2f deg: %s\n", alpha, msg);
            break
        end
    catch Esolve2
        fprintf("Failure detected, sweep stopped.\n")
        break;
    end
    Cls = [Cls, m.post.cl];
    Cds = [Cds, m.post.cd];
    Cms = [Cms, m.post.cm];
    alphas = [alphas alphas_neg(i)];
end

%% 
% Cleanup vectors
[alphas, I] = sort(alphas);
Cls = Cls(I);
Cds = Cds(I);
Cms = Cms(I);


%% Setup polar and slopes

if plotting == 1

    figure;
    t = tiledlayout(3,3);

    % Plot Airfoil
    ax1 = nexttile(1,[1 3]);
    plot(m.geom.xpoint(1,:),m.geom.xpoint(2,:))
    daspect(ax1,[1 1 1]);
    axis padded;
    title(airfoil,'interpreter','latex')

    % Plot Cl vs alpha
    ax2 = nexttile(4);
    plot(alphas, Cls)
    xlabel('$\alpha$ [deg]','interpreter','latex')
    ylabel('$C_l$','interpreter','latex')
    title('$C_l$ vs. $\alpha$','interpreter','latex')

    % Plot Cd vs alpha
    ax3 = nexttile(5);
    plot(alphas, Cds)
    xlabel('$\alpha$ [deg]','interpreter','latex')
    ylabel('$C_d$','interpreter','latex')
    title('$C_d$ vs. $\alpha$','interpreter','latex')

    % Plot Cm vs alpha
    ax4 = nexttile(6);
    plot(alphas, Cms)
    xlabel('$\alpha$ [deg]','interpreter','latex')
    ylabel('$C_m$','interpreter','latex')
    title('$C_m$ vs. $\alpha$','interpreter','latex')

    % Plot Cl vs Cd polar
    ax5 = nexttile(7);
    plot(Cds, Cls)
    xlabel('$C_d$','interpreter','latex')
    ylabel('$C_l$','interpreter','latex')
    title('Drag Polar','interpreter','latex')

    % Plot L'/D' vs Cl
    ax6 = nexttile(8);
    plot(Cls, Cls./Cds)
    xlabel('$C_l$','interpreter','latex')
    ylabel("$L'/D'$",'interpreter','latex')
    title("$L'/D'$ vs. $C_l$",'interpreter','latex')

    % Plot Cm vs Cl
    ax7 = nexttile(9);
    plot(Cls, Cms)
    xlabel('$C_l$','interpreter','latex')
    ylabel('$C_m$','interpreter','latex')
    title('$C_m$ vs. $C_l$','interpreter','latex')

end


%% Data output

% General
data.Cls = Cls;
data.Cds = Cds;
data.Cms = Cms;
data.alphas = alphas;

% Stall
data.Clmax = max(Cls);
data.alpha_max = alphas(Cls == data.Clmax);
data.Clmin = min(Cls);
data.alpha_min = alphas(Cls == data.Clmin);

% L'/D'
data.LDs = Cls./Cds;
data.LDmax = max(data.LDs);
data.alpha_LDmax = alphas(data.LDs == data.LDmax);


%% CDCL

% Find min Cd point
Cdmin = min(Cds);
Cl_Cdmin = Cls(Cds == Cdmin);

% Find stall points
Cd_Clmin = Cds(Cls == data.Clmin);
Cd_Clmax = Cds(Cls == data.Clmax);

% Make CDCL entry
CDCL = sprintf("%.4f %.4f   %.4f %.4f   %.4f %.4f   %.4f",...
    data.Clmin, Cd_Clmin, Cl_Cdmin, Cdmin, data.Clmax, Cd_Clmax);


%% CLAF 

% Find dCl/alpha [-/deg]
i_alpha0 = find(alphas == 0);
i_alpha5 = find(alphas == 3);
Cl_alpha = (Cls(i_alpha5) - Cls(i_alpha0))/(alphas(i_alpha5) - alphas(i_alpha0));

% Convert to dCl/alpha [-/rad]
Cl_alpha_rad = Cl_alpha*180/pid;

% Calculate CLAF factor (see AVL documentation)
CLAF = Cl_alpha_rad/(2*pi);


%% Save to .mat

name = strrep(string(airfoil), ' ', '_');
filename = sprintf('%s.mat',name);
save(filename,'-struct','data');
out = sprintf('\nAirfoil analysis saved to %s\n',filename);
fprintf(out)


end