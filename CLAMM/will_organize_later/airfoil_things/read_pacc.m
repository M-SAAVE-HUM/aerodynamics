function [polar] = read_pacc(filename)


    % Open file
    fid = fopen(filename,'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Read file line-by-line
    lines = {};
    while ~feof(fid)
        lines{end+1,1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % Collect numeric rows
    data = [];

    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue
        end

        nums = sscanf(line, '%f');

        % Expect at least 7 numeric columns for XFOIL polar data
        if numel(nums) >= 7
            data(end+1,:) = nums(1:7).'; %#ok<AGROW>
        end
    end

    if isempty(data)
        error('No polar data found in file: %s', filename);
    end

    % Extract columns
    alpha = data(:,1);
    cl    = data(:,2);
    cd    = data(:,3);
    cm    = data(:,5);

    % Sort by alpha
    [alpha_sorted, idx] = sort(alpha);

    polar.alpha = alpha_sorted;
    polar.cl    = cl(idx);
    polar.cd    = cd(idx);
    polar.cm    = cm(idx);


    %% CDCL

    % Find min Cd point
    Cdmin = min(cd);
    Cl_Cdmin = cl(cd == Cdmin);
    
    % Find stall points
    Cd_Clmin = cd(cl == min(cl));
    Cd_Clmax = cd(cl == max(cl));
    
    % Make CDCL entry
    CDCL = sprintf("%.4f %.4f   %.4f %.4f   %.4f %.4f   %.4f",...
    min(cl), Cd_Clmin, Cl_Cdmin, Cdmin, max(cl), Cd_Clmax);

    polar.CDCL = CDCL;


   %% CLAF

   % Find dCl/alpha [-/deg]
    i_alpha0 = find(alpha == 0,1);
    i_alpha5 = find(alpha == 3,1);
    Cl_alpha = (cl(i_alpha5) - cl(i_alpha0))/(alpha(i_alpha5) - alpha(i_alpha0));
    
    % Convert to dCl/alpha [-/rad]
    Cl_alpha_rad = Cl_alpha*180/pi;
    
    % Calculate CLAF factor (see AVL documentation)
    CLAF = Cl_alpha_rad/(2*pi);

    if CLAF < 1e-2
        CLAF = 0;
    end

    polar.CLAF = CLAF;


    %% Save to file

    name = strrep(string(filename), ' ', '_');
    name = strtrim(erase(name,'.txt'));
    filename = sprintf('%s.mat',name);
    save(filename,'-struct','polar');
    out = sprintf('\nAirfoil analysis saved to %s\n',filename);
    fprintf(out)

end

