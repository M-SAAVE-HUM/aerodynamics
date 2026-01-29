%% Lift Estimator

% inputs: aspect ratio (AR), speed (v), wing area (S), air density (rho),
% angle of attack (alpha), zero-lift angle of attack (a0) - in degrees

% outputs: lift coefficient (CL), lift force (N)

function[CL, L] = lift_estimator(AR, v, S, rho, alpha, a0)

    gamma = 1.4;                % specific heat ratio
    R = 287.185;                % specific gas constant (J/kgK)
    T = 288.15;                 % temperature (K)
    A = sqrt(gamma*R*T);        % speed of sound (m/s)
    M = v / A;                  % mach number
    alpha_rad = deg2rad(alpha); 
    a0_rad = deg2rad(a0);       
    a = (2*pi) / (1 + 2/AR); 

    % outputs
    CL = ( a / (sqrt(1 - M^2)) ) * (alpha_rad - a0_rad);
    L = 0.5 * rho * v^2 * S * CL; 

end