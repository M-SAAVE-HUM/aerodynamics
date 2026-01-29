% import RAE 2822 airfoil coordinates 
X = load('rae.txt');

% make/panel an airfoil out of these coordinates  
m = mfoil('coords', X, 'npanel',199);

% derotate: make the chord line horizontal 
m.geom_derotate
    
cls = zeros(1, 14);
cds = zeros(1, 14);
for i = 1:14
    % set operating conditions 
    m.setoper('alpha', i-4, 'Re', 6500000, 'Ma', 0.729);
    
    % run the solver, plot the results 
    m.solve 
    cls(i) = m.post.cl;
    cds(i) = m.post.cd;
end 

