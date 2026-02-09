%%  M-SAAVE CLAMM Sizing Constraint Plots
% Original code:
    % AE740 Final Project - Preliminary Sizing Plot
    % Kabir Khwaja, 4/17/25
% Modified:
    % Andrew Painter, 2/2/26

close all; clear; clc;

%% INPUTS

% Parameters
h = 256;            % altitude of Ann Arbor [m]
[~,~,~,rho] = atmoscoesa(h);        % air density [kg/m^3]
dTakeoff = 50;      % max takeoff distance [m]
dLand = 50;         % max landing distance [m]
Vcruise = 24;       % cruise speed [m/s]
Vturn = Vcruise;    % turning speed [m/s]
Vstall = 14;        % stall speed (desired) [m/s]
CLmax = 1.5;       % max. lift coefficient [-]
CD0 = 0.03;         % parasitic drag coefficient [-] (0.0132)
AR = 7.87;          % aspect ratio [-]
eCruise = 0.9;      % span efficiency factor at cruise [-]
eClimb = 0.85;      % span efficiency factor at climb [-]
G = 8;              % climb angle [degrees]
theta = 30;         % bank angle during turn maneuver [degrees]
g = 9.81;           % gravitational acceleration [m/s^2]


% Estimated design point #1 (based off Kabir's Python Aero Toolbox MDO results)
est1.T0 = 31;        % thrust (random guess) [N]
est1.m0 = 10.7;      % mass [kg]
est1.W0 = est1.m0*g; % weight [N]
est1.Sref0 = 0.8013; % wing area [m^2]
est1.WS0 = est1.W0/est1.Sref0; % wing loading [N/m^2]
est1.TW0 = est1.T0/est1.W0;    % thrust-to-weight ratio [N/N]
est1.color = 'r';    % dot color

% Estimated design point #2 (most W/S)
est2.T0 = 40;        % thrust (random guess) [N]
est2.m0 = 10.7;      % mass [kg]
est2.W0 = est2.m0*g; % weight [N]
est2.Sref0 = 0.6479; % wing area [m^2]
est2.WS0 = est2.W0/est2.Sref0; % wing loading [N/m^2]
est2.TW0 = est2.T0/est2.W0;    % thrust-to-weight ratio [N/N]
est2.color = 'g';    % dot color

% Estimated design point #3 (most W/S for least T/W)
est3.T0 = 26;        % thrust (random guess) [N]
est3.m0 = 10.7;      % mass [kg]
est3.W0 = est3.m0*g; % weight [N]
est3.Sref0 = 0.9719; % wing area [m^2]
est3.WS0 = est3.W0/est3.Sref0; % wing loading [N/m^2]
est3.TW0 = est3.T0/est3.W0;    % thrust-to-weight ratio [N/N]
est3.color = 'b';    % dot color




%%  Calculate constraints

const.WSvals = linspace(0, 280, 2801); % vector of wing loading values
const.WSstall = stall(rho, CLmax, Vstall);
const.WSlanding = landing(rho, dLand, CLmax);
const.TWtakeoff = takeoff(g, rho, dTakeoff, CLmax, const.WSvals);
const.TWclimb = climb(CLmax, CD0, AR, eClimb, G);
const.TWcruise = cruise(Vcruise, CD0, AR, eCruise, rho, const.WSvals);
const.TWceiling = ceiling(CD0, AR, eCruise);
const.TWmaneuver = maneuver(rho, theta, Vturn, AR, eCruise, CD0, const.WSvals);


%% T/W vs. W/S plot

fig1 = figure(1);
hold on
plot(const.WSvals, const.TWtakeoff,  'LineWidth',2, 'DisplayName','Takeoff')
plot(const.WSvals, const.TWmaneuver, 'LineWidth',2, 'DisplayName','Maneuver')
plot(const.WSvals, const.TWcruise,   'LineWidth',2, 'DisplayName','Cruise')
yline(const.TWclimb,   'LineWidth',2,'Color','#5c9c0e', 'DisplayName','Climb')
yline(const.TWceiling, 'LineWidth',2,'Color','#0e7d9c', 'DisplayName','Ceiling')
xline(const.WSstall,   'LineWidth',2,'Color','#7b0e9c', 'DisplayName','Stall Speed')
xline(const.WSlanding, 'LineWidth',2,'Color','#9c4c0e', 'DisplayName','Landing')

% Estimated design point
scatter(est1.WS0, est1.TW0, 80, est1.color, 'filled', 'DisplayName','Estimated Design Point #1')
scatter(est2.WS0, est2.TW0, 80, est2.color, 'filled', 'DisplayName','Estimated Design Point #2')
scatter(est3.WS0, est3.TW0, 80, est3.color, 'filled', 'DisplayName','Estimated Design Point #3')


% Fills
% constrained area to the right of the stall line/landing line
fill([const.WSstall 280 280 const.WSstall], [0 0 1 1], 'r', ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');
% fill([WSlanding 280 280 WSlanding], [0 0 1 1], 'r', ...
%     'FaceAlpha', 0.1, 'EdgeColor', 'none');

% area below each horizontal constraint line
fill([0 280 280 0], [0 0 const.TWclimb const.TWclimb], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');
fill([0 280 280 0], [0 0 const.TWceiling const.TWceiling], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');

% area below each curve-based constraint
fill([const.WSvals, fliplr(const.WSvals)], [const.TWtakeoff, zeros(size(const.WSvals))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');
WSvalsClip = const.WSvals(2:end);
maneuverConstraint = const.TWmaneuver(2:end);
cruiseConstraint = const.TWcruise(2:end);
zeroY = zeros(1, 2800);
fill([WSvalsClip, fliplr(WSvalsClip)], [maneuverConstraint, zeroY], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');
fill([WSvalsClip, fliplr(WSvalsClip)], [cruiseConstraint, zeroY], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off');

% Axes
xlabel('W/S [$N/m^2$]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('T/W', 'Interpreter', 'latex', 'FontSize', 16)
legend('Location','northeast');
xlim([0 280])
ylim([0 1])
hold off

% saveas(fig1, 'TW_WS_sizing_plot.eps', 'epsc');


%% T/S plots

[TSfig1] = TS_plot(est1,const);
[TSfig2] = TS_plot(est2,const);
[TSfig3] = TS_plot(est3,const);

% saveas(TSfig1, 'T_S_sizing_plot.eps', 'epsc');
% saveas(TSfig2, 'T_S_sizing_plot.eps', 'epsc');
% saveas(TSfig3, 'T_S_sizing_plot.eps', 'epsc');


%% T/S plot function

function [fig] = TS_plot(est,const)

    West = est.W0;

    % Convert to S and T
    Svals = (const.WSvals/West).^(-1);
    Slanding = (const.WSlanding/West).^(-1);
    Sstall = (const.WSstall/West).^(-1);
    Ttakeoff = const.TWtakeoff.*West;
    Tclimb = const.TWclimb.*West;
    Tcruise = const.TWcruise.*West;
    Tceiling = const.TWceiling.*West;
    Tmaneuver = const.TWmaneuver.*West;
    
    fig = figure;
    hold on
    plot(Svals, Ttakeoff,  'LineWidth',2, 'DisplayName','Takeoff')
    plot(Svals, Tmaneuver, 'LineWidth',2, 'DisplayName','Maneuver')
    plot(Svals, Tcruise,   'LineWidth',2, 'DisplayName','Cruise')
    yline(Tclimb,   'LineWidth',2, 'Color','#5c9c0e', 'DisplayName','Climb')
    yline(Tceiling, 'LineWidth',2, 'Color','#0e7d9c', 'DisplayName','Ceiling')
    xline(Sstall,   'LineWidth',2, 'Color','#7b0e9c', 'DisplayName','Stall')
    xline(Slanding, 'LineWidth',2, 'Color','#9c4c0e', 'DisplayName','Landing')
    
    % Estimated design point
    scatter(est.Sref0, est.T0, 80, 'r', 'filled', 'DisplayName','Estimated Design Point')
    
    
    xlabel('$S_{ref}$ [$m^2$]','Interpreter','Latex','FontSize', 16)
    ylabel('$T$ [$N$]','Interpreter','Latex','FontSize', 16)
    legend('Location','northeast');
    xlim([0 6])
    ylim([0 400])
    str = sprintf('Estimated Mass = %2.2f kg',est.m0);
    title(str,'FontSize', 16)
    hold off
    
end


%% Stall Function
% INPUTS:  rho - air density [kg/m^3]
%          CLmax - estimated max lift coefficient [-]
%          Vstall - desired stall speed [m/s]
% OUTPUTS: WSstall - wing loading for stall constraint [N/m^2]

function [WSstall] = stall(rho, CLmax, Vstall)
    WSstall = 0.5*rho*(Vstall.^2)*CLmax;
end

%% Takeoff Function
% INPUTS:  g - gravitational acceleration [m/s^2]
%          rho - air density [kg/m^3]
%          dTakeoff - max. takeoff distance [m]
%          CLmax - max. lift coefficient [-]
%          WS - vector of wing loading values [N/m^2]
% OUTPUTS: TWtakeoff - min. thrust to weight ratio required for takeoff [-]

function [TWtakeoff] = takeoff(g, rho, dTakeoff, CLmax, WS)
    k = 1.1; % takeoff speed to stall speed ratio
    SF = 1.5; % safety factor to consider drag/other deficiencies
    TWtakeoff = (SF.*WS.*k^2) ./ (g*rho*dTakeoff*CLmax);
end

%% Landing Function
% INPUTS:  rho - air density [kg/m^3]
%          dLand - max. landing distance [m]
%          CLmax - max. lift coefficient [-]
% OUTPUTS: WSlanding - wing loading for landing constraint [N/m^2]
% TODO - figure out how to convert this properly
function [WSlanding] = landing(rho, dLand, CLmax)
    % sigma = 0.794;  % landing/sea-level density for hot day at 5000 ft
    dObst = 0;  % 50 ft tall obstacle avoidance distance in meters [m]
    sigma = 1;
    % WSlanding = 47.88*(dLand - dObst)*(sigma)*CLmax/80;
    % WSlanding = (dLand*1.5 + 1000)/(25*rho*CLmax); % raymer formula converted to m
    CLmax = 2.5;
    WSlanding = dLand*rho*CLmax/(0.592*2);
end



%% Climb Function
% INPUTS:  CLmax - max lift coefficient [-]
%          CD0 - parasitic drag [-]
%          G - climb angle [degrees]
%          AR - aspect ratio [-]
%          e - span efficiency factor at climb [-]
% OUTPUTS: TWclimb - min. thrust to weight ratio needed for climb [-]

function [TWclimb] = climb(CLmax, CD0, AR, e, G)
    ks = 1.2; % ratio of climb speed to stall speed
    CLclimb = CLmax/(ks^2); % CL at climb
    G = tand(G); % climb gradient percent
    TWclimb = ((CD0 + CLclimb^2/(pi*AR*e)) / CLclimb) + G; 
end

%% Cruise Function
% INPUTS:  Vcruise - cruise speed [m/s]
%          CD0 - parasitic drag coefficient [-]
%          AR - aspect ratio [-]
%          e - span efficiency at cruise [-]
%          rho - air density [kg/m^3]
%          WS - vector of wing loading values [N/m^2]
% OUTPUTS: TWcruise - thrust to weight ratio needed for cruise [-]

function [TWcruise] = cruise(Vcruise, CD0, AR, e, rho, WS)
    q = 0.5*rho*(Vcruise^2); % dynamic pressure [Pa]
    TWcruise = ((q*CD0)./WS) + (WS./(q*pi*AR*e)); 
end

%% Ceiling Function
% INPUTS:  CD0 - parasitic drag coefficient [-]
%          AR - aspect ratio [-]
%          e - span efficiency at cruise [-]
% OUTPUTS: TWceiling - thrust to weight ratio at ceiling

function [TWceiling] = ceiling(CD0, AR, e)
    G = 0.001; % small climb gradient correction
    K = 1/(pi*AR*e);
    TWceiling = G + 2*sqrt(CD0*K);
end

%% Maneuver Function
% INPUTS:  rho - air density [kg/m^3]
%          theta - bank angle [degrees]
%          Vturn - turning speed [m/s]
%          AR - aspect ratio [-]
%          e - span efficiency at cruise [-]
%          CD0 - parasitic drag coefficient [-]
%          WS - vector of wing loading values [N/m^2]
% OUTPUTS: TWmaneuver - thrust to weight needed for banked turn [-]

function [TWmaneuver] = maneuver(rho, theta, Vturn, AR, e, CD0, WS)
    n = 1/cosd(theta); % load factor
    q = 0.5*rho*Vturn^2; % dynamic pressure
    TWmaneuver = (q*CD0)./WS + WS.*(n^2 / (q*pi*AR*e));
end
