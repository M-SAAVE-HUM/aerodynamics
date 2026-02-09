function [] = gen_mass(ac,rho,name,folder)

mass = ac.mass;
xcg = ac.xcg;
ycg = ac.ycg;
zcg = ac.zcg;
output = string(sprintf(['# Mass properties file for AVL\n' ...
    '# Only total mass and Xcg defined\n' ...
    'Lunit = 1 m\n' ...
    'Munit = 1 kg\n' ...
    'Tunit = 1 s\n'...
    'g = 9.81\n'...
    'rho = %.6f\n'...
    '# mass (kg), x, y, z, Ixx, Iyy, Izz, Ixy, Ixz, Iyz\n'...
    '%.6f %.6f %.6f %.6f 0 0 0 0 0 0\n'],rho,mass,xcg,ycg,zcg));


%% Output file

if ~isfolder(folder)
    mkdir(folder);
end

filename = folder + "/" + name + ".mass";
fid = fopen(filename,'w');
fprintf(fid,"%s", output);
fclose(fid);

out = sprintf('\tMass file saved to %s\n',filename);
fprintf(out)


end