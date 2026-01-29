%% Zero Lift Calculator!!!
% use this to find zero lift angle of attack for a given airfoil!!!
% runs an airfoil in mfoil and enforces a CL of 0

% load airfoil
%X = load('ch10.txt');
m = mfoil('naca', '6412', 'npanel', 199);

% set coefficient of lift to zero  
m.setoper('cl', 0.0);

% run the solver, plot the results 
% alpha = zero-lift angle of attack
m.solve 
%m.oper
