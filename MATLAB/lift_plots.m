%% lift plots

% gonk parameters
AR = 6.5; % aspect ratio 
S = 1.5; % wing area (m^2)
a0 = -6.2957; % zero lift aoa (degrees)
rho = 1.225; % air density (kg/m^3)

% set up matrix
K = 61;
velocities = linspace(0,24,K);
alphas = linspace(-5,10,K);
cl_speed_alpha = zeros(K, K);
lift_speed_alpha = zeros(K, K);

% loop to add values to matrix
for i = 1:K
    for j = 1:K
        [CL, L] = lift_estimator(AR, velocities(i), S, rho, alphas(j), a0);
        cl_speed_alpha(i,j) = CL;
        lift_speed_alpha(i,j) = L;
    end
end

zero_alpha_index = find(alphas==0);
lift_speed = zeros(1, K);
lift_alpha = zeros(1, K);
for i = 1:K
    lift_speed(i) = lift_speed_alpha(i, zero_alpha_index);
    lift_alpha(i) = lift_speed_alpha(K,i);
end

%% plots
t = tiledlayout(2,3);

% lift vs speed
nexttile
plot(velocities, lift_speed, '-r', 'Linewidth', 1.5)
xlabel('Velocity (m/s)')
ylabel('Lift (N)')
title('Lift vs Velocity (at Angle of Attack = 0')
grid on;

% Span across two rows and columns
nexttile([2 2])
% surface plot
[X,Y] = meshgrid(velocities, alphas);
surf(X, Y, lift_speed_alpha);
xlabel('Velocity (m/s)')
ylabel('Angle of Attack (degrees)')
zlabel('Lift (N)')
title('Gonk Lift vs Operating Velocities and Angles of Attack')

% Last tile
nexttile
plot(alphas, lift_alpha, '-b', 'Linewidth', 1.5)
xlabel('Angle of Attack (degrees)')
ylabel('Lift (N)')
title('Lift vs Angle of Attack (at Cruise Speed = 24 m/s)')
grid on;