%% Drag Estimator
% made with love by Kabir <3
% fixes to make: add in landing gear/prop drag estimates (see Hoerner)

% notes:
% numbers fitted for gonk aircraft, form factors from Raymer
% total drag coefficient, CD = CD_parasitic + CD_induced = CD0 + CDi

% inputs: air density (kg/m^3), speed (m/s), Reynold's number, Oswald
% efficiency factor, aspect ratio, wing area (m^2), mass (kg)

% outputs: coeff. parasitic drag, coeff. induced drag, coeff. total drag, 
% and total drag force

function[CD0, CDi, CD, D] = drag_estimator(rho, u, Re, e, AR, S, m)

    g = 9.81; % gravity (m/s^2)
    W = m * g; % weight (N)
    CL = W/(0.5 * rho * u^2 * S); % lift coefficient
    Cf = 0.0576 / (Re)^0.2; % skin friction coefficient
    K = 1/(pi*AR*e); % induced drag factor
    Q = 2; % interference (fudge) factor

    % wing
    FF_wing = 1 + (0.6/0.3)*(0.12) + 100*(0.12)^4; % form factor wing
    S_wetwing = 2 * S; % wing wetted area

    % tail
    FF_tail = 1 + (0.6/0.3)*(0.12) + 100*(0.12)^4;
    S_wettail = 2 * (0.309 + 0.134);

    % fuselage
    w = 0.283; % approx width of fuselage (m)
    h = 0.260; % approx height of fuselage (m)
    d = 0.275; % approx diameter of fuselage (m)
    l = 1.826; % approx length fuselage (nose to taper end, in m)
    ffr = l/d; % fuselage fineness ratio
    FF_fuselage = 0.9 + (5 / (ffr^1.5)) + (ffr / 400);
    front_area = w * h; % frontal area
    sides = 4 * l * d; % total area of sides
    S_wetfuse = front_area + sides;    

    % outputs
    CD0 = (Cf*Q/S) * ( (FF_wing * S_wetwing) + (FF_tail * S_wettail) ...
        + (FF_fuselage * S_wetfuse)); % parasitic drag

    % ^ compare this to CD0 from NX surface area
    CDi = K*(CL)^2; % induced drag
    CD = CD0 + CDi; % total craft drag
    D = 0.5*rho*S*CD*u^2; % drag

end