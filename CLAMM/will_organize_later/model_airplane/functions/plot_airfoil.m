function [outputArg1,outputArg2] = plot_airfoil(filenames)
% Plots important variables for the given airfoil


for i = 1:length(filenames)
    af{i} = load(filenames(i));
end

af{2}


%% Setup polar and slopes

figure;
t = tiledlayout(3,3);

% Plot Airfoil
ax1 = nexttile(1,[1 3]);
hold on
for i = 1:length(af)
    label = sprintf("%s, Re = %d, M = %.2f",string(af{i}.name), af{i}.Re, af{i}.Mach);
    plot(0,0,'displayName',label)

end
daspect(ax1,[1 1 1]);
legend();
axis padded;
hold off

% Plot Cl vs alpha
ax2 = nexttile(4);
hold on
for i = 1:length(af)
    plot(af{i}.alpha, af{i}.cl)
end
xlabel('$\alpha$ [deg]','interpreter','latex')
ylabel('$C_l$','interpreter','latex')
title('$C_l$ vs. $\alpha$','interpreter','latex')
hold off

% Plot Cd vs alpha
hold on
ax3 = nexttile(5);
for i = 1:length(af)
    af{1}.cd
    af{1}.alpha
    plot(af{i}.alpha, af{i}.cd)
end
xlabel('$\alpha$ [deg]','interpreter','latex')
ylabel('$C_d$','interpreter','latex')
title('$C_d$ vs. $\alpha$','interpreter','latex')
hold off

% Plot Cm vs alpha
hold on
ax4 = nexttile(6);
for i = 1:length(af)
    plot(af{i}.alpha, af{i}.cm)
end
xlabel('$\alpha$ [deg]','interpreter','latex')
ylabel('$C_m$','interpreter','latex')
title('$C_m$ vs. $\alpha$','interpreter','latex')
hold off

% Plot Cl vs Cd polar
ax5 = nexttile(7);
hold on
for i = 1:length(af)
    plot(af{i}.cd, af{i}.cl)
end
xlabel('$C_d$','interpreter','latex')
ylabel('$C_l$','interpreter','latex')
title('Drag Polar','interpreter','latex')
hold off

% Plot L'/D' vs Cl
ax6 = nexttile(8);
hold on
for i = 1:length(af)
    plot(af{i}.cl, af{i}.cl./af{i}.cd)
end
xlabel('$C_l$','interpreter','latex')
ylabel("$L'/D'$",'interpreter','latex')
title("$L'/D'$ vs. $C_l$",'interpreter','latex')
hold off

% Plot Cm vs Cl
ax7 = nexttile(9);
hold on
for i = 1:length(af)
    plot(af{i}.cl, af{i}.cm)
end
xlabel('$C_l$','interpreter','latex')
ylabel('$C_m$','interpreter','latex')
title('$C_m$ vs. $C_l$','interpreter','latex')
hold off




%% Helper

function [x, y] = readAirfoilDat(filename)
% readAirfoilDat  Read x,y coordinates from an airfoil .dat file

    data = readmatrix(filename, ...
        'FileType','text', ...
        'Delimiter',{' ','\t'}, ...
        'ConsecutiveDelimitersRule','join');

    data = data(~any(isnan(data),2), :); % Remove any non-numerics

    x = data(:,1);
    y = data(:,2);
end

end