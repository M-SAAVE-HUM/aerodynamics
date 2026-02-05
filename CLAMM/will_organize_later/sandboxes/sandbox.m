clear
clc
close all



af = load("2412_xfoil.mat");

alphas = af.alpha;
Cls = af.cl;
Cds = af.cd;
Cms = af.cm;


plotting = 1;
if plotting == 1

    figure;
    t = tiledlayout(3,3);

    % Plot Airfoil
    ax1 = nexttile(1,[1 3]);
    plot(0,0)
    daspect(ax1,[1 1 1]);
    axis padded;
    title('shutup','interpreter','latex')

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


%%

wing.Sref = 0.8013;       % Reference planform area [m^2]
wing.AR = 7.87;           % Aspect Ratio [-]
wing.taper1 = 1.0;        % Taper ratio 1 (c_root/c_mid-station) [-]
wing.taper2 = 0.5;        % Taper ratio 2 (c_mid-station/c_tip) [-]
wing.btaper = 0.5;        % Half-span fraction of mid-station
wing.dihedral1 = 0.0;     % Dihedral of root-mid-span section [deg]
wing.dihedral2 = 0.0;     % Dihedral of mid-span-tip section [deg]
wing.incidence = 0.0;     % Angle of incidence of wing relative to fuselage centerline [deg]
wing.xLE = 0.0;           % Root leading edge x-location
wing.yLE = 0.0;           % Root leading edge y-location
wing.zLE = 0.0;           % Root leading edge z-location

wing.aroot = 'naca2412.dat';
wing.CLAF_root = af.CLAF;
wing.CDCL_root = af.CDCL;

wing.akink = 'naca2412.dat';
wing.CLAF_kink = af.CLAF;
wing.CDCL_kink = af.CDCL;

wing.atip = 'naca2412.dat';
wing.CLAF_tip = af.CLAF;
wing.CDCL_tip = af.CDCL;

name = 'test';
template = 'template_3.avl';
folder = '\wings\';

gen_wing(wing,name,template,folder)



