%% airfoil comparison

% Sweep through NACA max cambers (first digit)
alphas = linspace(0, 10, 11);
num_airfoils = 4;
cls_camber = zeros(num_airfoils, length(alphas));
cds_camber = zeros(num_airfoils, length(alphas));

num = 2412; % starting airfoil number

for i = 1:num_airfoils
    % Convert number to NACA string
    number_str = sprintf('%04d', num);

    % Create airfoil object
    m = mfoil('naca', number_str, 'npanel', 199);

    % Sweep through alpha
    for j = 1:length(alphas)
        m.setoper('alpha', alphas(j), 'Re', 800000);
        m.solve;
        cls_camber(i, j) = m.post.cl;
        cds_camber(i, j) = m.post.cd;
    end

    % Increment to next airfoil
    num = num + 2000;
end
 
% sweep through max camber locations
cls_locs = zeros(num_airfoils, length(alphas));
cds_locs = zeros(num_airfoils, length(alphas));

num = 6212; % starting airfoil number

for i = 1:num_airfoils
    % Convert number to NACA string
    number_str = sprintf('%04d', num);

    % Create airfoil object
    m = mfoil('naca', number_str, 'npanel', 199);

    % Sweep through alpha
    for j = 1:length(alphas)
        m.setoper('alpha', alphas(j), 'Re', 800000);
        m.solve;
        cls_locs(i, j) = m.post.cl;
        cds_locs(i, j) = m.post.cd;
    end

    % Increment to next airfoil
    num = num + 200;
end

% Plot results
figure(2);
hold on;
for k = 1:num_airfoils
    plot(alphas, cls_camber(k,:), 'DisplayName', sprintf('NACA %04d', 1000*k + 412), 'LineWidth', 2);
end
xlabel('Angle of Attack (°)');
ylabel('C_l');
legend show;

figure(3);
hold on;
for k = 1:num_airfoils
    plot(alphas, cls_locs(k,:), 'DisplayName', sprintf('NACA %04d', 6012 + 200*k), 'LineWidth', 2);
end
xlabel('Angle of Attack (°)');
ylabel('C_l');
legend show;

figure(4);
hold on;
for k = 1:num_airfoils
    plot(cds_camber(k,:), cls_camber(k,:), 'DisplayName', sprintf('NACA %04d', 6012 + 200*k), 'LineWidth', 2);
    plot(cds_locs(k,:), cls_locs(k,:), 'DisplayName', sprintf('NACA %04d', 6012 + 200*k), 'LineWidth', 2);
end
xlabel('C_d');
ylabel('C_l');
legend show;

% sweep through thicknesses