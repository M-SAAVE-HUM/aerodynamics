%% Mission Profile Plotter
% plots of altitude and location for mission duration

% phases of flight
% 1. takeoff

syms v
CD_stall

dc = @(v) .5.*CDstall.*rho.*(v^2).*Sref;   
Vcruise = fzero(@(v) Tcurve(v) - dc(v), 1);