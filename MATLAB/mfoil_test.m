% start with a 4-digit NACA airfoil 
m = mfoil('naca','6412', 'npanel',199);

% set angle of attack and Reynolds number  
m.setoper('alpha', 10, 'Re', 10^6);

% run the solver, plot the results 
m.solve 