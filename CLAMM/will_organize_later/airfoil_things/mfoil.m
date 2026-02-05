%-------------------------------------------------------------------------------
% mfoil.m: class and methods for subsonic airfoil analysis (v 2025-04-24)
%
% Copyright (C) 2025 Krzysztof J. Fidkowski
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%-------------------------------------------------------------------------------

classdef mfoil < handle
  
%===============================================================================
  properties

    version = '2025-04-28';    % mfoil version
    geom  = struct_geom();     % geometry
    foil  = struct_panel();    % airfoil panels
    wake  = struct_panel();    % wake panels
    oper  = struct_oper();     % operating conditions
    isol  = struct_isol();     % inviscid solution variables
    vsol  = struct_vsol();     % viscous solution variables
    glob  = struct_glob();     % global system variables
    post  = struct_post();     % post-processing quantities
    param = struct_param();    % parameters
    
  end
  
  
%===============================================================================
  methods

    %-------------------------------------------------------
    % Constructor function
    function M = mfoil(varargin)
% Commented out for compatibility with earlier Matlab releases
% $$$       arguments
% $$$         args.coords (:,:) {mustBeNumeric} = []      % coordinate matrix
% $$$         args.naca   (1,:) char            = '0012'  % NACA digits
% $$$         args.npanel (1,1) {mustBeInteger} = 199     % number of panels
% $$$       end
      
      p = inputParser;
      addParameter(p, 'naca', '0012');
      addParameter(p, 'coords', []);
      addParameter(p, 'npanel', 199);
      parse(p, varargin{:});
      coords = p.Results.coords;
      naca = p.Results.naca;
      npanel = p.Results.npanel;
      
      if (~isempty(coords))
        set_coords(M, coords);
      else
        naca_points(M, naca);
      end
      make_panels(M, npanel);
    end 
    
    %-------------------------------------------------------
    % Set operating conditions
    function setoper(M, varargin)
% Commented out for compatibility with earlier Matlab releases
% $$$       arguments
% $$$         M          (1,1) mfoil
% $$$         args.alpha (1,1) {mustBeNumeric}  % angle of attack (deg)
% $$$         args.cl    (1,1) {mustBeNumeric}  % target lift coefficient
% $$$         args.Re    (1,1) {mustBeNumeric}  % Reynolds number
% $$$         args.Ma    (1,1) {mustBeNumeric}  % Mach number
% $$$         args.xftl  (1,1) {mustBeNumeric}  % forced transition x/c, lower
% $$$         args.xftu  (1,1) {mustBeNumeric}  % forced transition x/c, upper
% $$$         args.visc  (1,1) {mustBeNumericOrLogical} % true -> viscous
% $$$       end

      p = inputParser;
      addParameter(p, 'alpha', nan);
      addParameter(p, 'cl', nan);
      addParameter(p, 'Re', nan);
      addParameter(p, 'Ma', nan);
      addParameter(p, 'visc', nan);
      addParameter(p, 'xftl', nan);
      addParameter(p, 'xftu', nan);
      parse(p, varargin{:});
      args = struct;
      if (~isnan(p.Results.alpha)), args.alpha = p.Results.alpha; end;
      if (~isnan(p.Results.cl)),    args.cl    = p.Results.cl;    end;
      if (~isnan(p.Results.Re)),    args.Re    = p.Results.Re;    end;
      if (~isnan(p.Results.Ma)),    args.Ma    = p.Results.Ma;    end;
      if (~isnan(p.Results.xftl)),  args.xftl  = p.Results.xftl;  end;
      if (~isnan(p.Results.xftu)),  args.xftu  = p.Results.xftu;  end;
      if (~isnan(p.Results.visc)),  args.visc  = p.Results.visc;  end;
      
      
      if isfield(args, 'alpha') 
        M.oper.alpha = args.alpha; M.oper.givencl = false;
      end
      if isfield(args, 'cl')
        M.oper.cltgt = args.cl; M.oper.givencl = true;
      end
      if isfield(args, 'Re')
        M.oper.Re = args.Re; M.oper.viscous = true;
      end
      if isfield(args, 'Ma'), M.oper.Ma = args.Ma; end
      if isfield(args, 'xftl')
        M.oper.xft(1) = args.xftl;
      end
      if isfield(args, 'xftu')
        M.oper.xft(2) = args.xftu;
      end
      if isfield(args, 'visc')
        if (M.oper.viscous ~= args.visc), clear_solution(M); end;
        M.oper.viscous = args.visc; 
      end
    end
    
    %-------------------------------------------------------
    % Solve current point
    function solve(M)
      if (M.oper.viscous)
        solve_viscous(M);
        if (M.param.doplot), close all; mplot_results(M); end
      else 
        solve_inviscid(M);
        if (M.param.doplot), close all; mplot_results(M); end
      end
    end
    
    %-------------------------------------------------------
    % Geometry functions
    function geom_change(M, args) % change airfoil coords
      arguments
        M           (1,1) mfoil
        args.coords (:,:) {mustBeNumeric} = []      % coordinate matrix
        args.naca   (1,:) char            = '0012'  % NACA digits
      end
      if (~isempty(args.coords))
        set_coords(M, args.coords);
      else
        naca_points(M, args.naca);
      end
      make_panels(M, M.foil.N-1);
    end
    function geom_flap(M, xzhinge, eta)
      mgeom_flap(M, xzhinge, eta); % add a flap
    end
    function geom_addcamber(M, zcamb)
      mgeom_addcamber(M, zcamb); % increment camber
    end
    function geom_derotate(M)
      mgeom_derotate(M); % derotate: make chordline horizontal
    end
    function geom_panel(M, npanel)
      make_panels(M, npanel); % panel the airfoil
    end
    
    %-------------------------------------------------------
    % Plotting functions
    function plot_results(M)
      mplot_results(M);  % summary results plot
    end   
    function plot_panels(M)
      mplot_panels(M);  % airfoil/wake paneling plot
    end   
    function plot_airfoil(M)
      mplot_airfoil(M);  % airfoil geometry plot
    end   
    function plot_cpplus(M)
      mplot_cpplus(M);  % cp and output plot
    end   
    function plot_boundary_layer(M)
      mplot_boundary_layer(M);  % boundary layer plot
    end   
    function plot_stream(M, xrange, npoint)
      mplot_stream(M, xrange, npoint);  % streamfunction
    end
    function plot_velocity(M, xrange, npoint)
      mplot_velocity(M, xrange, npoint); % velocity_vectors
    end
    function plot_distributions(M)
      mplot_distributions(M); % viscous distributions
    end
  
    
    %-------------------------------------------------------
    % Verify derivatives by finite differences
    function ping(M)
      %system_test(M);
      ping_test(M);
    end
    
  end
  
end


%==================   STRUCTURES   ===================

%-------------------------------------------------------------------------------
function S = struct_geom()
  S.chord = 1;        % chord length
  S.wakelen = 1.0;    % wake extent length, in chords
  S.npoint = 0;       % number of geometry representation points
  S.name = 'noname';  % airfoil name, e.g. NACA XXXX
  S.xpoint = [];      % point coordinates, [2 x npoint]
  S.xref = [0.25; 0]; % moment reference point
end

%-------------------------------------------------------------------------------
function S = struct_panel()
  S.N = 0;  % number of nodes
  S.x = []; % node coordinates, [2 x N]
  S.s = []; % arclength values at nodes
  S.t = []; % dx/ds, dy/ds tangents at nodes
end

%-------------------------------------------------------------------------------
function S = struct_oper()
  S.Vinf = 1;          % velocity magnitude
  S.alpha = 0;         % angle of attack, in degrees
  S.rho = 1;           % density
  S.cltgt = 0;         % lift coefficient target
  S.givencl = false;   % true if cl is given instead of alpha
  S.initbl = true;     % true to initialize the boundary layer
  S.viscous = false;   % true to do viscous
  S.redowake = false;  % true to rebuild wake after alpha changes
  S.Re = 1e5;          % viscous Reynolds number
  S.Ma = 0;            % Mach number
  S.forcet = [0,0,0];  % forced transition flags (internal use only)
  S.xft = [1,1,1];     % forced transition x/c values (lower, upper, wake)
  S.xift = [0,0,0];    % forced transition xi values (internal use only)
end

%-------------------------------------------------------------------------------
function S = struct_isol()
  S.AIC = [];          % aero influence coeff matrix
  S.gamref = [];       % 0,90-deg alpha vortex strengths at airfoil nodes
  S.gam = [];          % vortex strengths at airfoil nodes (for current alpha)
  S.sstag = 0.;        % s location of stagnation point
  S.sstag_g = [0,0];   % lin of sstag w.r.t. adjacent gammas
  S.sstag_ue = [0,0];  % lin of sstag w.r.t. adjacent ue values
  S.Istag = [0,0];     % node indices before/after stagnation
  S.sgnue = [];        % +/- 1 on upper/lower surface nodes
  S.xi = [];           % distance from the stagnation at all points
  S.uewi = [];         % inviscid edge velocity in wake
  S.uewiref = [];      % 0,90-deg alpha inviscid ue solutions on wake
end

%-------------------------------------------------------------------------------
function S = struct_vsol()
  S.th  = [];       % theta = momentum thickness [Nsys]
  S.ds = [];        % delta star = displacement thickness [Nsys]
  S.Is = {};        % 3 cell arrays of surface indices
  S.wgap = [];      % wake gap over wake points
  S.ue_m = [];      % linearization of ue w.r.t. mass (all nodes)
  S.sigma_m = [];   % d(source)/d(mass) matrix
  S.ue_sigma = [];  % d(ue)/d(source) matrix
  S.turb = [];      % flag over nodes indicating if turbulent (1) or lam (0) 
  S.xt = 0.;        % transition location (xi) on current surface under consideration
  S.Xt = [0,0;0,0]; % transition xi/x for lower and upper surfaces 
end

%-------------------------------------------------------------------------------
function S = struct_glob()
  S.Nsys = 0;      % number of equations and states
  S.U = [];        % primary states (th,ds,sa,ue) [4 x Nsys]
  S.dU = [];       % primary state update
  S.dalpha = 0;    % angle of attack update
  S.conv = true;   % converged flag
  S.R = [];        % residuals [3*Nsys x 1]
  S.R_U = [];      % residual Jacobian w.r.t. primary states
  S.R_x = [];      % residual Jacobian w.r.t. xi (s-values) [3*Nsys x Nsys]
end

%-------------------------------------------------------------------------------
function S = struct_post()
  S.cp = [];       % cp distribution
  S.cpi = [];      % inviscid cp distribution
  S.cl = 0;        % lift coefficient
  S.cl_ue = [];    % linearization of cl w.r.t. ue [N, airfoil only]
  S.cl_alpha = 0;  % linearization of cl w.r.t. alpha
  S.cm = 0;        % moment coefficient
  S.cdpi = 0;      % near-field pressure drag coefficient
  S.cd = 0;        % total drag coefficient
  S.cdf = 0;       % skin friction drag coefficient
  S.cdp = 0;       % pressure drag coefficient

  % distributions
  S.th = [];  % theta = momentum thickness distribution
  S.ds = [];  % delta* = displacement thickness distribution
  S.sa = [];  % amplification factor/shear lag coeff distribution
  S.ue = [];  % edge velocity (compressible) distribution
  S.uei = []; % inviscid edge velocity (compressible) distribution
  S.cf = [];  % skin friction distribution
  S.Ret = []; % Re_theta distribution
  S.Hk = [];  % kinematic shape parameter distribution
end

%-------------------------------------------------------------------------------
function S = struct_param()

  S.verb   = 1;     % printing verbosity level (higher -> more verbose)
  S.rtol   = 1e-10; % residual tolerance for Newton
  S.niglob = 50;    % maximum number of global iterations
  S.doplot = true;  % true to plot results after each solution
  S.axplot = [];    % plotting axes (for more control of where plots go)
  S.legpos = [.7,.15]; % cp plot legend position (x,y) coords
  S.vlegpos = [.3,.88]; % viscous legend position (x,y) coords
  
  % viscous parameters
  S.ncrit  = 9.0;   % critical amplification factor    
  S.Cuq    = 1.0;   % scales the uq term in the shear lag equation
  S.Dlr    = 0.9;   % wall/wake dissipation length ratio
  S.SlagK  = 5.6;   % shear lag constant
  
  % initial Ctau after transition
  S.CtauC  = 1.8;   % Ctau constant
  S.CtauE  = 3.3;   % Ctau exponent
  
  % G-beta locus: G = GA*sqrt(1+GB*beta) + GC/(H*Rt*sqrt(cf/2))
  S.GA     = 6.7;   % G-beta A constant
  S.GB     = 0.75;  % G-beta B constant
  S.GC     = 18.0;  % G-beta C constant
  
  % operating conditions and thermodynamics
  S.Minf   = 0.;    % freestream Mach number
  S.Vinf   = 0.;    % freestream speed
  S.muinf  = 0.;    % freestream dynamic viscosity
  S.mu0    = 0.;    % stagnation dynamic viscosity
  S.rho0   = 1.;    % stagnation density
  S.H0     = 0.;    % stagnation enthalpy
  S.Tsrat  = 0.35;  % Sutherland Ts/Tref
  S.gam    = 1.4;   % ratio of specific heats
  S.KTb    = 1.;    % Karman-Tsien beta = sqrt(1-Minf^2)
  S.KTl    = 0.;    % Karman-Tsien lambda = Minf^2/(1+KTb)^2
  S.cps    = 0.;    % sonic cp
  
  % station information
  S.simi   = false; % true at a similarity station
  S.turb   = false; % true at a turbulent station
  S.wake   = false; % true at a wake station
end


%% ============ INPUT AND OUTPUT ==============


%-------------------------------------------------------------------------------
function vprint(param, verb, varargin)
% prints according to a verbosity flag (higher -> more verbose)
  if (verb <= param.verb), fprintf(1,varargin{:}); end
end


%% ============ PLOTTING AND POST-PROCESSING  ==============


%-------------------------------------------------------------------------------
function mplot_panels(M)
% plots the airfoil and wake panels
% INPUT
%   M : mfoil structure
% OUTPUT
%   plot of panels
  ax = M.param.axplot; if isempty(ax), ax = gca; end;
  plot(ax, M.foil.x(1,:), M.foil.x(2,:), 'bo-', 'linewidth', 2); hold on;
  if (size(M.wake.x,2) > 0)
    plot(ax, M.wake.x(1,:), M.wake.x(2,:), 'ro-', 'linewidth', 2);
  end
end


%-------------------------------------------------------------------------------
function mplot_stream(M, xrange, npoint)
% plots the streamfunction contours
% INPUT
%   M : mfoil structure
%   xrange = [xmin, xmax, zmin, zmax] window range
%   npoint = number of 1d plotting points
% OUTPUT
%   streamfunction plot
  assert(~isempty(M.isol.gam), 'No inviscid solution');
  if isempty(xrange), xrange = [-.05, 1.05, -.1, .1]; end
  sx = linspace(xrange(1), xrange(2), npoint);
  sz = linspace(xrange(3), xrange(4), npoint);
  
  [Z,X] = meshgrid(sz,sx); P = X;
  
  N = M.foil.N;         % number of airfoil points
  Nw = M.wake.N;        % number of wake points
  Vinf = M.oper.Vinf;   % freestream speed
  alpha = M.oper.alpha; % angle of attack [deg]
  G = M.isol.gam;       % gamma vector

  [~,~,~,tcp,tdp] = TE_info(M.foil.x); % trailing-edge info
  
  % calculate source terms if viscous
  if (M.oper.viscous)
    ue = M.glob.U(4,:)'; ds = M.glob.U(2,:)';
    sigv = M.vsol.sigma_m*(ue.*ds);
  end
  
  % loop over plotting points
  for ip = 1:size(X,1)
    for jp = 1:size(X,2)
      
      xi = [X(ip,jp), Z(ip,jp)];
      
      Psi = 0;
      
      % panel influence
      for j=1:N-1
        [aij, bij] = panel_linvortex_stream(M.foil.x(:,[j,j+1]), xi);
        Psi = Psi + aij*G(j) + bij*G(j+1);
      end
      
      % TE source influence
      a = panel_constsource_stream(M.foil.x(:,[N,1]), xi);
      Psi = Psi + a*0.5*tcp*(G(N)-G(1));

      % TE vortex influence
      [a,b] = panel_linvortex_stream(M.foil.x(:,[N,1]), xi);
      Psi = Psi + (a+b)*0.5*(-tdp)*(G(N)-G(1));
      
      % viscous source influence
      if (M.oper.viscous)
        for i = 1:N-1
          Psi = Psi + sigv(i)*panel_constsource_stream(M.foil.x(:,[i,i+1]), xi);
        end
        for i = 1:Nw-1
          Psi = Psi + sigv(N-1+i)*panel_constsource_stream(M.wake.x(:,[i,i+1]), xi);
        end
      end
      
      % add freestream
      P(ip,jp) = Psi + Vinf*cosd(alpha)*Z(ip,jp) - Vinf*sind(alpha)*X(ip,jp);
    end
  end
  
  figure; clf;
  contourf(X,Z,P,51); hold on;
  plot_panels(M);
  set('fontsize', 14);
  axis equal; axis(xrange);
  
end


%-------------------------------------------------------------------------------
function mplot_velocity(M, xrange, npoint)
% makes a quiver plot of velocity vectors in the domain
% INPUT
%   xrange : axes range [xmin, xmax, zmin, zmax]
%   npoint : number of points (1d)
% OUTPUT
%   figure with velocity vectors
  
  N = M.foil.N; Nw = M.wake.N;  % number of points on the airfoil/wake
  if isempty(xrange), xrange = [-.05, 1.05, -.1, .1]; end
  sx = linspace(xrange(1), xrange(2), npoint);
  sz = linspace(xrange(3), xrange(4), npoint);
  
  [Z,X] = meshgrid(sz,sx); U = X; V = Z;
  % calculate source terms if viscous
  if (M.oper.viscous)
    ue = M.glob.U(4,:)'; ds = M.glob.U(2,:)';
    sv = M.vsol.sigma_m*(ue.*ds);
  end
  for ip = 1:npoint
    for jp = 1:npoint
      % global coordinate
      xij = [X(ip,jp), Z(ip,jp)];
      % first the inviscid velocity
      v = inviscid_velocity(M.foil.x, M.isol.gam, M.oper.Vinf, M.oper.alpha, xij);
      if (M.oper.viscous)
        % next, the viscous source contribution
        for i = 1:N-1
          sigma = sv(i);
          v = v + sigma*panel_constsource_velocity(M.foil.x(:,[i,i+1]), xij, []);
        end
        for i = 1:Nw-1
          sigma = sv(N-1+i); %dm/ds;
          v = v + sigma*panel_constsource_velocity(M.wake.x(:,[i,i+1]), xij, []);
        end
      end
      U(ip,jp) = v(1); V(ip,jp) = v(2);
    end
  end
  figure(1); clf;
  quiver(X,Z,U,V); hold on;
  plot_panels(M);
  set(gca, 'fontsize', 14);
  axis equal
  axis(xrange);
  
end


%-------------------------------------------------------------------------------
function mplot_cpplus(M)
% makes a cp plot with outputs printed
% INPUT
%   M : mfoil structure
% OUTPUT
%   cp plot on current axes

  ax = M.param.axplot; if isempty(ax), ax = gca; end;
  %fg = M.param.fgplot; if isempty(fg), fg = gcf; end;
  chord = M.geom.chord;
  xz = [M.foil.x, M.wake.x]; x = xz(1,:); N = M.foil.N;
  xrange = [-.1,1.4]*chord;
  if (M.oper.viscous > 0)
    ctype = {'r', 'b', 'k'};
    for is = 1:3
      Is = M.vsol.Is{is};
      plot(ax, x(Is), M.post.cp(Is), [ctype{is},'-'], 'linewidth', 2); hold(ax,'on');
      plot(ax, x(Is), M.post.cpi(Is), [ctype{is},':'], 'linewidth', 2); 
    end
  else
    plot(ax, x, M.post.cp, 'b-', 'linewidth', 2); hold(ax,'on');
  end
  if (M.oper.Ma > 0) && (M.param.cps > (min(M.post.cp)-.2))
    plot(ax, [xrange(1), chord], M.param.cps*[1,1], 'k--', 'linewidth', 2);
    text(ax, 0.8*chord, M.param.cps-0.1, 'sonic $c_p$', 'interpreter', 'latex', 'fontsize', 18);
  end

  set(ax, 'fontsize', 16); box(ax,'off');
  ylabel(ax, '$c_p$', 'interpreter', 'latex', 'position', [-0.18, 0, -1], 'fontsize', 24);
  V = axis(ax); V(1:2) = xrange; V(4) = 1; axis(ax,V);
  set(ax, 'xtick', linspace(0,1.4*chord, 15));
  set(ax, 'ydir', 'reverse');   % invert the cp axis
  set(ax, 'color', 'none');
  set(ax, 'xaxislocation', 'origin')
  set(ax, 'xticklabel', [])
  
  % output text box
  name = sprintf('%s', M.geom.name); 
  if isempty(name), name = 'Untitled'; end

  str = {sprintf('\\underline{%s}', name), ...
         sprintf('Ma = $%.4f$', M.oper.Ma), ...
         sprintf('$\\alpha = %.4f^{\\circ}$', M.oper.alpha), ...
         sprintf('$c_{\\ell} = %.4f$', M.post.cl), ...
         sprintf('$c_{m} = %.4f$', M.post.cm), ...
         sprintf('$c_{d} = %.4f$', M.post.cd)};
  lp = M.param.legpos;
  xloc = ax.XLim*[1-lp(1);lp(1)]; yloc = ax.YLim*[1-lp(2); lp(2)];
  h = text(ax, xloc, yloc, str, 'interpreter', 'latex', 'fontsize', 16);
  if (M.oper.viscous > 0)
    str = {sprintf('Re = %.1e', M.oper.Re), ...
           sprintf('$c_{df} = %.5f$', M.post.cdf), ...
           sprintf('$c_{dp} = %.5f$', M.post.cdp)};
    if (M.oper.xft(1)<1)
      str = {str{:}, sprintf('tripl $x/c$ = %.3f', M.oper.xft(1))};
    end
    if (M.oper.xft(2)<1)
      str = {str{:}, sprintf('tripu $x/c$ = %.3f', M.oper.xft(2))};
    end
    lp = M.param.vlegpos;
    xloc = ax.XLim*[1-lp(1);lp(1)]; yloc = ax.YLim*[1-lp(2); lp(2)];
    h = text(ax, xloc, yloc, str, 'interpreter', 'latex', 'fontsize', 16);
  end

end


%-------------------------------------------------------------------------------
function mplot_airfoil(M)
% makes an airfoil plot
% INPUT
%   M : mfoil structure
% OUTPUT
%   airfoil plot on current axes

  ax = M.param.axplot; if isempty(ax), ax = gca; end;
  chord = M.geom.chord;
  xz = [M.foil.x, M.wake.x]; x = xz(1,:); N = M.foil.N;
  xrange = [-.1,1.4]*chord;
  plot(ax, xz(1,:), xz(2,:), 'k-', 'linewidth', 1); hold(ax,'on');
  P = pbaspect(ax); AR = P(2);
  V = axis(ax); V(1:2) = xrange; 
  dy = V(4)-V(3); dx = xrange(2)-xrange(1); V(3:4) = V(3:4)*(dx*AR)/dy;
  axis(ax,V); box(ax, 'off'); axis(ax,'off');
  
end

%-------------------------------------------------------------------------------
function mplot_boundary_layer(M)
% makes a plot of the boundary layer
% INPUT
%   M : mfoil structure
% OUTPUT
%   boundary layer plot on current axes

  if (M.oper.viscous <= 0), return; end
  ax = M.param.axplot; if isempty(ax), ax = gca; end;  
  xz = [M.foil.x, M.wake.x]; x = xz(1,:); N = M.foil.N;
  ds = M.post.ds; % displacement thickness
  rl = 0.5*(1+(ds(1)-ds(N))/ds(N+1)); ru = 1-rl;
  t = [M.foil.t, M.wake.t]; % tangents
  n = [-t(2,:); t(1,:)]; n = n./vecnorm(n); % outward normals
  xzd = xz + n.*ds; % airfoil + delta*
  ctype = {'r', 'b', 'k'};
  for i = 1:4
    is = i;
    if (is==3), xzd = xz + n.*ds*ru; end
    if (is==4), xzd = xz - n.*ds*rl; is = 3; end
    Is = M.vsol.Is{is};
    plot(ax, xzd(1,Is), xzd(2,Is), [ctype{is},'-'], 'linewidth', 2);
    hold(ax, 'on');
  end
  
end


%-------------------------------------------------------------------------------
function mplot_results(M)
% makes a summary results plot with cp, airfoil, BL delta, outputs
% INPUT
%   M : mfoil structure
% OUTPUT
%   summary results plot as a new figure  
    
  assert(~isempty(M.post.cp), 'no cp for results plot');
  figure;
  
  % pressure coefficient
  SS = get(groot, 'screensize'); % to make an well-proportioned figure
  npx = ceil(min(SS(3:4))/2); vs = [npx, ceil(npx*0.8)];
  V = get(gcf, 'position'); V(3:4) = vs; set(gcf, 'position', V);
  
  % cp plot
  subplot('position', [0.1 0.3, 0.85, 0.65]);
  mplot_cpplus(M);
  
  % airfoil plot
  subplot('position', [0.1 0.05, 0.85, 0.22]);
  mplot_airfoil(M);
    
  % BL thickness
  mplot_boundary_layer(M);
  
end


%-------------------------------------------------------------------------------
function calc_force(M)
% calculates force and moment coefficients
% INPUT
%   M : mfoil structure with solution (inviscid or viscous)
% OUTPUT
%   M.post values are filled in
% DETAILS
%   lift/moment are computed from a panel pressure integration
%   the cp distribution is stored as well
  
  chord = M.geom.chord; xref = M.geom.xref; % chord and ref moment point 
  Vinf = M.param.Vinf; rho = M.oper.rho; alpha = M.oper.alpha;
  qinf = 0.5*rho*Vinf^2; % dynamic pressure
  N = M.foil.N; % number of points on the airfoil
  
  % calculate the pressure coefficient at each node
  if (M.oper.viscous), ue = M.glob.U(4,:); else, ue = get_ueinv(M)'; end
  [cp, cp_ue] = get_cp(ue, M.param); M.post.cp = cp;
  M.post.cpi = get_cp(get_ueinv(M)',M.param); % inviscid cp
  
  % lift, moment, near-field pressure cd coefficients by cp integration  
  cl = 0; cl_ue = zeros(1,N); cl_alpha = 0; cm = 0; cdpi = 0;  
  for i0 = 2:N+1
    i = i0; ip=i-1; if (i0==N+1), i=1; ip=N; end
    x1 = M.foil.x(:,ip); x2 = M.foil.x(:,i); % panel points
    dxv = x2-x1; dx1 = x1-xref; dx2 = x2-xref;
    dx1nds = dxv(1)*dx1(1)+dxv(2)*dx1(2); % (x1-xref) cross n*ds
    dx2nds = dxv(1)*dx2(1)+dxv(2)*dx2(2); % (x2-xref) cross n*ds
    dx = -dxv(1)*cosd(alpha) - dxv(2)*sind(alpha); % minus from CW node ordering
    dz =  dxv(2)*cosd(alpha) - dxv(1)*sind(alpha); % for drag
    cp1 = cp(ip); cp2 = cp(i); cpbar = 0.5*(cp(ip)+cp(i)); % average cp on the panel
    cl = cl + dx*cpbar;
    I = [ip,i]; cl_ue(I) = cl_ue(I) + dx*0.5*cp_ue(I);
    cl_alpha = cl_alpha + cpbar*(sind(alpha)*dxv(1) - cosd(alpha)*dxv(2))*pi/180;
    cm = cm + cp1*dx1nds/3 + cp1*dx2nds/6 + cp2*dx1nds/6 + cp2*dx2nds/3;
    cdpi = cdpi + dz*cpbar;
  end
  cl = cl/chord; cm = cm/chord^2; cdpi = cdpi/chord;
  M.post.cl = cl; M.post.cl_ue = cl_ue; M.post.cl_alpha = cl_alpha;
  M.post.cm = cm; M.post.cdpi = cdpi;
  
  % viscous contributions
  cd = 0; cdf = 0;
  if (M.oper.viscous)
    
    % Squire-Young relation for total drag (exrapolates theta from end of wake)
    iw = M.vsol.Is{3}(end); % station at the end of the wake
    U = M.glob.U(:,iw); H = U(2)/U(1); ue = get_uk(U(4), M.param); % state
    cd = 2.0*U(1)*(ue/Vinf)^((5+H)/2.);
    
    % skin friction drag
    Df = 0.;
    for is = 1:2
      Is = M.vsol.Is{is}; % surface point indices
      param = build_param(M, is); % get parameter structure
      param = station_param(M, param, Is(1));
      cf1 = 0; %get_cf(M.glob.U(:,Is(1)), param); % first cf value
      ue1 = 0; %get_uk(M.glob.U(4,Is(1)), param);
      rho1 = rho;
      x1 = M.isol.xstag;
      for i = 1:length(Is) % loop over points
        param = station_param(M, param, Is(i));
        cf2 = get_cf(M.glob.U(:,Is(i)), param); % get cf value
        ue2 = get_uk(M.glob.U(4,Is(i)), param);
        rho2 = get_rho(M.glob.U(:,Is(i)), param);
        x2 = M.foil.x(:,Is(i)); dxv = x2 - x1;
        dx = dxv(1)*cosd(alpha) + dxv(2)*sind(alpha);
        Df = Df + 0.25*(rho1*cf1*ue1^2 + rho2*cf2*ue2^2)*dx;
        cf1 = cf2; ue1 = ue2; x1 = x2; rho1 = rho2;
      end
    end
    cdf = Df/(qinf*chord);    
  end
  % store results
  M.post.cd = cd; M.post.cdf = cdf; M.post.cdp = cd-cdf;

  % print out current values
  fmt = '  alpha=%.2fdeg, cl=%.6f, cm=%.6f, cdpi=%.6f, cd=%.6f, cdf=%.6f, cdp=%.6f\n';
  vprint(M.param,1, fmt, M.oper.alpha, M.post.cl, M.post.cm, M.post.cdpi, ...
         M.post.cd, M.post.cdf, M.post.cdp);
    
end


%-------------------------------------------------------------------------------
function get_distributions(M)
% computes various distributions (quantities at nodes) and stores them in M.post
% INPUT
%   M  : mfoil class with a valid solution in M.glob.U
% OUTPUT
%   M.post : distribution quantities calculated
% DETAILS
%   Relevant for viscous solutions
  
  assert(~isempty(M.glob.U), 'no global solution');
  
  % quantities already in the global state
  M.post.th = M.glob.U(1,:); % theta
  M.post.ds = M.glob.U(2,:); % delta*
  M.post.sa = M.glob.U(3,:); % amp or ctau
  M.post.ue = get_uk(M.glob.U(4,:), M.param); % compressible edge velocity 
  M.post.uei = get_ueinv(M); % compressible inviscid edge velocity
  
  % derived viscous quantities
  N = M.glob.Nsys; cf = zeros(N,1); Ret = zeros(N,1); Hk = zeros(N,1);
  for is = 1:3   % loop over surfaces
    Is = M.vsol.Is{is}; % surface point indices
    param = build_param(M, is); % get parameter structure
    for i = 1:length(Is)  % loop over points
      j = Is(i); Uj = M.glob.U(:,j);
      param = station_param(M, param, j);
      uk = get_uk(Uj(4), param); % corrected edge speed
      cfloc = get_cf(Uj, param); % local skin friction coefficient
      cf(j) = cfloc * uk^2/param.Vinf^2; % free-stream-based cf
      Ret(j) = get_Ret(Uj, param); % Re_theta
      Hk(j) = get_Hk(Uj, param); % kinematic shape factor
    end
  end
  M.post.cf = cf; M.post.Ret = Ret; M.post.Hk = Hk;
  
end


%-------------------------------------------------------------------------------
function plot_quantity(M, q, qname)
% plots a quantity q over lower/upper/wake surfaces
% INPUT
%   M     : mfoil class with valid airfoil/wake points
%   q     : vector of values to plot, on all points (wake too if present)
%   qname : name of quantity, for axis labeling
% OUTPUT
%   figure showing q versus x
  
  figure;
  xy = [M.foil.x, M.wake.x];  % xy = M.isol.xi; % uncomment to plot vs xi
  ctype = {'r-', 'b-', 'k-'};
  if (isempty(M.vsol.Is)), plot(xy(1,:), q, [ctype{1},'o'], 'linewidth', 2);
  else
    sleg = {'lower', 'upper', 'wake'};
    for is = 1:3
      Is = M.vsol.Is{is};
      plot(xy(1,Is), q(Is), [ctype{is},'o'], 'linewidth', 2, 'DisplayName', sleg{is});
      hold on;
    end
  end
  set(gca, 'fontsize', 16); grid on;
  xlabel('x = distance along chord', 'fontsize', 18);
  ylabel(qname, 'fontsize', 18);
  if (~isempty(M.vsol.Is)),h = legend('Location', 'SouthEast'); set(h, 'fontsize', 18);  end
  %set(gcf, 'Position', [2800, 2100-300*fi, 1000, 300]);

end


%-------------------------------------------------------------------------------
function plot_mass(M)
% plots viscous/mass/source quantities
% INPUT
%   M     : mfoil class with valid solution
% OUTPUT
%   Plots of mass flow at each node, and source on each panel
 
  N = M.foil.N; Nw = M.wake.N;  % number of points on the airfoil/wake
  xpan = [0.5*(M.foil.x(1,1:N-1)+M.foil.x(1,2:N)), 0.5*(M.wake.x(1,1:Nw-1)+M.wake.x(1,2:Nw))];
  ue = M.glob.U(4,:)'; ds = M.glob.U(2,:)'; 
  m = ue.*ds; % mass flow at nodes
  sigma = M.vsol.sigma_m*m; % source on panels

  plot_quantity(M, m, 'mass flow');
  figure;
  sleg = {'lower', 'upper', 'wake'};
  ctype = {'b-', 'r-', 'k-'};
  is = 1; Is = 1:(N-1); 
  plot(xpan(Is), sigma(Is), [ctype{is},'o'], 'linewidth', 2, 'DisplayName', sleg{is});
  hold on;
  is = 3; Is = N:(N+Nw-2); 
  plot(xpan(Is), sigma(Is), [ctype{is},'o'], 'linewidth', 2, 'DisplayName', sleg{is});
  set(gca, 'fontsize', 16); grid on;
  xlabel('\xi = distance from stagnation point', 'fontsize', 18);
  ylabel('source', 'fontsize', 18);

end


%-------------------------------------------------------------------------------
function mplot_distributions(M)
% plots viscous distributions
% INPUT
%   M  : mfoil class with solution
% OUTPUT
%   figures of viscous distributions

  close all;
  get_distributions(M);
  plot_quantity(M, M.post.ue, 'u_e = edge velocity');
  plot_quantity(M, M.post.uei, 'inviscid edge velocity');
  plot_quantity(M, M.post.sa, 'amplification or c_{\tau}^{1/2}');
  plot_quantity(M, M.post.Hk, 'H_k = kinematic shape parameter');
  plot_quantity(M, M.post.ds, '\delta^* = displacement thickness');
  plot_quantity(M, M.post.th, '\theta = momentum thickness');
  plot_quantity(M, M.post.cf, 'c_f = skin friction coefficient');
  plot_quantity(M, M.post.Ret, 'Re_{\theta} = theta Reynolds number');
  plot_mass(M);
  
end


%% ============ INVISCID FUNCTIONS ==============


%-------------------------------------------------------------------------------
function clear_solution(M)
% clears inviscid/viscous solutions by re-initializing structures
% INPUT
%   M : mfoil structure
% OUTPUT
%   M : mfoil structure without inviscid or viscous solution
% DETAILS

  M.isol = struct_isol();
  M.vsol = struct_vsol();
  M.glob = struct_glob();
  M.post = struct_post();
  M.wake.N = 0;
  M.wake.x = [];
  M.wake.s = [];
  M.wake.t = [];

end


%-------------------------------------------------------------------------------
function solve_inviscid(M)
% solves the inviscid system, rebuilds 0,90deg solutions
% INPUT
%   M : mfoil structure
% OUTPUT
%   inviscid vorticity distribution is computed
% DETAILS
%   Uses the angle of attack in M.oper.gamma
%   Also initializes thermo variables for normalization

  assert(M.foil.N>0, 'No panels');
  M.oper.viscous = false;
  init_thermo(M);
  M.isol.sgnue = ones(1,M.foil.N); % do not distinguish sign of ue if inviscid
  build_gamma(M, M.oper.alpha);
  if (M.oper.givencl), cltrim_inviscid(M); end
  calc_force(M);
  M.glob.conv = true; % no coupled system ... convergence is guaranteed

end


%-------------------------------------------------------------------------------
function cltrim_inviscid(M)
% trims inviscid solution to prescribed target cl, using alpha
% INPUT
%   M : mfoil structure
% OUTPUT
%   inviscid vorticity distribution is computed for a given cl
% DETAILS
%   Iterates using cl_alpha computed in post-processing
%   Accounts for cl_ue in total derivative

  for i = 1:15 % trim iterations
    alpha = M.oper.alpha; % current angle of attack
    calc_force(M); % calculate current cl and linearization
    R = M.post.cl - M.oper.cltgt;
    if (norm(R) < 1e-10), break; end
    sc = [-sind(alpha);cosd(alpha)]*pi/180;
    cl_a = M.post.cl_alpha + M.post.cl_ue*(M.isol.gamref*sc); % total deriv
    dalpha = -R/cl_a;
    M.oper.alpha = alpha + min(max(dalpha,-2), 2);
  end
  if (i>=15), vprint(M.param,1, '** inviscid cl trim not converged **'); end
  M.isol.gam = M.isol.gamref*[cosd(M.oper.alpha); sind(M.oper.alpha)];
  
end


%-------------------------------------------------------------------------------
function [ueinv] = get_ueinv(M)
% computes invicid tangential velocity at every node
% INPUT
%   M : mfoil structure
% OUTPUT
%   ueinv : inviscid velocity at airfoil and wake (if exists) points
% DETAILS
%   The airfoil velocity is computed directly from gamma
%   The tangential velocity is measured + in the streamwise direction

  assert(~isempty(M.isol.gam), 'No inviscid solution');
  alpha = M.oper.alpha; cs = [cosd(alpha); sind(alpha)];
  uea = M.isol.sgnue'.*(M.isol.gamref*cs);  % airfoil
  if (M.oper.viscous) && (M.wake.N > 0)
    uew = M.isol.uewiref*cs; % wake
    uew(1) = uea(end); % ensures continuity of upper surface and wake ue
  else, uew = [];
  end
  ueinv = [uea; uew]; % airfoil/wake edge velocity
  
end


%-------------------------------------------------------------------------------
function [ueinvref] = get_ueinvref(M)
% computes 0,90deg inviscid tangential velocities at every node
% INPUT
%   M : mfoil structure
% OUTPUT
%   ueinvref : 0,90 inviscid tangential velocity at all points (N+Nw)x2
% DETAILS
%   Uses gamref for the airfoil, uewiref for the wake (if exists)
  
  assert(~isempty(M.isol.gam), 'No inviscid solution');
  uearef = M.isol.sgnue'.*M.isol.gamref; % airfoil
  if (M.oper.viscous) && (M.wake.N > 0)
    uewref = M.isol.uewiref; % wake
    uewref(1,:) = uearef(end,:); % continuity of upper surface and wake
  else, uewref = []; 
  end
  ueinvref = [uearef; uewref];
  
end


%-------------------------------------------------------------------------------
function build_gamma(M, alpha)
% builds and solves the inviscid linear system for alpha=0,90,input
% INPUT
%   M     : mfoil structure
%   alpha : angle of attack (degrees)
% OUTPUT
%   M.isol.gamref : 0,90deg vorticity distributions at each node (Nx2)
%   M.isol.gam    : gamma for the particular input angle, alpha
%   M.isol.AIC    : aerodynamic influence coefficient matrix, filled in
% DETAILS
%   Uses streamfunction approach: constant psi at each node
%   Continuous linear vorticity distribution on the airfoil panels
%   Enforces the Kutta condition at the TE
%   Accounts for TE gap through const source/vorticity panels
%   Accounts for sharp TE through gamma extrapolation

  N = M.foil.N;         % number of points  
  A = zeros(N+1,N+1);  % influence matrix
  rhs = zeros(N+1,2);  % right-hand sides for 0,90
  [~,hTE,~,tcp,tdp] = TE_info(M.foil.x); % trailing-edge info
  nogap = (abs(hTE) < 1e-10*M.geom.chord); % indicates no TE gap
  
  vprint(M.param,1, '\n <<< Solving inviscid problem >>> \n');
  
  % Build influence matrix and rhs
  for i=1:N             % loop over nodes
    xi = M.foil.x(:,i); % coord of node i
    for j=1:N-1         % loop over panels
      [aij, bij] = panel_linvortex_stream(M.foil.x(:,[j,j+1]), xi);
      A(i,j  ) = A(i,j  ) + aij;
      A(i,j+1) = A(i,j+1) + bij;
      A(i,N+1) = -1; % last unknown = streamfunction value on surf
    end
    % right-hand sides
    rhs(i,:) = [-xi(2), xi(1)];
    % TE source influence
    a = panel_constsource_stream(M.foil.x(:,[N,1]), xi);
    A(i,1) = A(i,1) - a*(0.5*tcp);
    A(i,N) = A(i,N) + a*(0.5*tcp);
    % TE vortex panel
    [a, b] = panel_linvortex_stream(M.foil.x(:,[N,1]), xi);
    A(i,1) = A(i,1) - (a+b)*(-0.5*tdp);
    A(i,N) = A(i,N) + (a+b)*(-0.5*tdp);
  end
  
  % special Nth equation (extrapolation of gamma differences) if no gap
  if (nogap), A(N,:) = 0; A(N,[1,2,3,N-2,N-1,N]) = [1,-2,1,-1,2,-1]; end

  % Kutta condition
  A(N+1,1) = 1;
  A(N+1,N) = 1;

  % Solve system for unknown vortex strengths
  M.isol.AIC = A;
  g = M.isol.AIC\rhs;
  M.isol.gamref = g(1:end-1,:); % last value is surf streamfunction
  M.isol.gam = M.isol.gamref(:,1)*cosd(alpha) + M.isol.gamref(:,2)*sind(alpha);

end


%-------------------------------------------------------------------------------
function varargout = inviscid_velocity(X, G, Vinf, alpha, x)
% returns inviscid velocity at x due to gamma (G) on panels X, and Vinf
% INPUT
%   X     : coordinates of N panel nodes (N-1 panels) (Nx2)
%   G     : vector of gamma values at each airfoil node (Nx1)
%   Vinf  : freestream speed magnitude
%   alpha : angle of attack (degrees)
%   x     : location of point at which velocity vector is desired  
% OUTPUT
%   V    : velocity at the desired point (2x1)
%   V_G  : (optional) linearization of V w.r.t. G, (2xN)
% DETAILS
%   Uses linear vortex panels on the airfoil
%   Accounts for TE const source/vortex panel
%   Includes the freestream contribution
  
  N = size(X,2);   % number of points  
  V = zeros(2,1);  % velocity
  dolin = (nargout > 1); % linearization requested
  if (dolin), V_G = zeros(2,N); end
  [~,~,~,tcp,tdp] = TE_info(X); % trailing-edge info
  % assume x is not a midpoint of a panel (can check for this)
  for j=1:(N-1)               % loop over panels
    [a, b] = panel_linvortex_velocity(X(:,[j,j+1]), x, [], false);
    V = V + a*G(j) + b*G(j+1);
    if (dolin), V_G(:,j) = V_G(:,j) + a; V_G(:,j+1) = V_G(:,j+1) + b; end
  end
  % TE source influence
  a = panel_constsource_velocity(X(:,[N,1]), x, []);
  f1 = a*(-0.5*tcp); f2 = a*0.5*tcp;
  V = V + f1*G(1) + f2*G(N);
  if (dolin), V_G(:,1) = V_G(:,1) + f1; V_G(:,N) = V_G(:,N) + f2; end
  % TE vortex influence
  [a,b] = panel_linvortex_velocity(X(:,[N,1]), x, [], false);
  f1 = (a+b)*(0.5*tdp); f2 = (a+b)*(-0.5*tdp);
  V = V + f1*G(1) + f2*G(N);
  if (dolin), V_G(:,1) = V_G(:,1) + f1; V_G(:,N) = V_G(:,N) + f2; end
  % freestream influence
  V = V + Vinf*[cosd(alpha); sind(alpha)];
  if (nargout > 0), varargout{1} = V; end
  if (nargout > 1), varargout{2} = V_G; end
  
end


%-------------------------------------------------------------------------------
function build_wake(M)
% builds wake panels from the inviscid solution
% INPUT
%   M     : mfoil class with a valid inviscid solution (gam)
% OUTPUT
%   M.wake.N  : Nw, the number of wake points
%   M.wake.x  : coordinates of the wake points (2xNw)
%   M.wake.s  : s-values of wake points (continuation of airfoil) (1xNw)
%   M.wake.t  : tangent vectors at wake points (2xNw)
% DETAILS
%   Constructs the wake path through repeated calls to inviscid_velocity
%   Uses a predictor-corrector method
%   Point spacing is geometric; prescribed wake length and number of points
  
  assert(~isempty(M.isol.gam), 'No inviscid solution');
  N = M.foil.N;  % number of points on the airfoil
  Vinf = M.oper.Vinf;    % freestream speed
  Nw = ceil(N/10 + 10*M.geom.wakelen); % number of points on wake
  S = M.foil.s;  % airfoil S values
  ds1 = 0.5*(S(2)-S(1) + S(N)-S(N-1)); % first nominal wake panel size
  sv = space_geom(ds1, M.geom.wakelen*M.geom.chord, Nw); % geometrically-spaced points
  xyw = zeros(2,Nw); tw = xyw; % arrays of x,y points and tangents on wake
  xy1 = M.foil.x(:,1); xyN = M.foil.x(:,N); % airfoil TE points
  xyte = 0.5*(xy1 + xyN); % TE midpoint
  n = xyN-xy1; t = [n(2); -n(1)]; % normal and tangent
  assert(t(1) > 0, 'Wrong wake direction; ensure airfoil points are CCW');
  xyw(:,1) = xyte + 1e-5*t*M.geom.chord; % first wake point, just behind TE
  sw = S(N) + sv; % s-values on wake, measured as continuation of the airfoil
  
  % loop over rest of wake
  for i = 1:Nw-1
    v1 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw(:,i));
    v1 = v1/norm(v1); tw(:,i) = v1; % normalized
    xyw(:,i+1) = xyw(:,i) + (sv(i+1)-sv(i))*v1; % forward Euler (predictor) step
    v2 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw(:,i+1));
    v2 = v2/norm(v2); tw(:,i+1) = v2; % normalized
    xyw(:,i+1) = xyw(:,i) + (sv(i+1)-sv(i))*0.5*(v1+v2); % corrector step
  end

  % determine inviscid ue in the wake, and 0,90deg ref ue too
  uewi = zeros(Nw,1); uewiref = zeros(Nw,2);
  for i = 1:Nw
    v = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw(:,i));
    uewi(i) = v'*tw(:,i);
    v = inviscid_velocity(M.foil.x, M.isol.gamref(:,1), Vinf, 0, xyw(:,i));
    uewiref(i,1) = v'*tw(:,i);
    v = inviscid_velocity(M.foil.x, M.isol.gamref(:,2), Vinf, 90, xyw(:,i));
    uewiref(i,2) = v'*tw(:,i);
  end
  
  % set values
  M.wake.N = Nw;
  M.wake.x = xyw;
  M.wake.s = sw;
  M.wake.t = tw;
  M.isol.uewi = uewi;
  M.isol.uewiref = uewiref;
  
end


%-------------------------------------------------------------------------------
function stagpoint_find(M)
% finds the LE stagnation point on the airfoil (using inviscid solution)
% INPUTS
%   M  : mfoil class with inviscid solution, gam
% OUTPUTS
%   M.isol.sstag   : scalar containing s value of stagnation point
%   M.isol.sstag_g : linearization of sstag w.r.t gamma (1xN)
%   M.isol.Istag   : [i,i+1] node indices before/after stagnation (1x2)
%   M.isol.sgnue   : sign conversion from CW to tangential velocity (1xN)
%   M.isol.xi      : distance from stagnation point at each node (1xN)

  assert(~isempty(M.isol.gam), 'No inviscid solution');
  N = M.foil.N;  % number of points on the airfoil
  J = find(M.isol.gam>0);  assert(~isempty(J), 'no stagnation point');
  I = [J(1)-1, J(1)]; G = M.isol.gam(I); S = M.foil.s(I);
  M.isol.Istag = I;  % indices of neighboring gammas
  den = (G(2)-G(1)); w1 = G(2)/den; w2 = -G(1)/den;
  sst = w1*S(1) + w2*S(2);  % s location
  M.isol.sstag = sst;
  M.isol.xstag = M.foil.x(:,I)*[w1;w2];  % x location
  st_g1 = G(2)*(S(1)-S(2))/(den*den);
  M.isol.sstag_g = [st_g1, -st_g1];
  sgnue = -1*ones(1,N); sgnue(J) = 1; % upper/lower surface sign
  M.isol.sgnue = sgnue;
  M.isol.xi = [abs(M.foil.s-M.isol.sstag), M.wake.s-M.isol.sstag];

end


%-------------------------------------------------------------------------------
function rebuild_isol(M)
% rebuilds inviscid solution, after an angle of attack change
% INPUT
%   M     : mfoil class with inviscid reference solution and angle of attack
% OUTPUT
%   M.isol.gam : correct combination of reference gammas
%   New stagnation point location if inviscid
%   New wake and source influence matrix if viscous

  assert(~isempty(M.isol.gam), 'No inviscid solution to rebuild');
  vprint(M.param,2, '\n  Rebuilding the inviscid solution.\n');
  alpha = M.oper.alpha;
  M.isol.gam = M.isol.gamref(:,1)*cosd(alpha) + M.isol.gamref(:,2)*sind(alpha);
  if (~M.oper.viscous)
    % viscous stagnation point movement is handled separately
    stagpoint_find(M);
  elseif (M.oper.redowake)
    build_wake(M);
    identify_surfaces(M);
    calc_ue_m(M); % rebuild matrices due to changed wake geometry
  end
  
end


%% ============ PANELING  ==============


%-------------------------------------------------------------------------------
function make_panels(M, npanel)
% places panels on the current airfoil, as described by M.geom.xpoint
% INPUT
%   M      : mfoil class
%   npanel : number of panels
% OUTPUT
%   M.foil.N : number of panel points
%   M.foil.x : coordinates of panel nodes (2xN)
%   M.foil.s : arclength values at nodes (1xN)
%   M.foil.t : tangent vectors, not normalized, dx/ds, dy/ds (2xN)
% DETAILS
%   Uses curvature-based point distribution on a spline of the points

  clear_solution(M); % clear any existing solution
  Ufac = 2;  % uniformity factor (higher, more uniform paneling)
  TEfac = 0.1; % Trailing-edge factor (higher, more TE resolution)
  [M.foil.x, M.foil.s, M.foil.t] = spline_curvature(M.geom.xpoint, npanel+1, Ufac, TEfac);
  M.foil.N = size(M.foil.x,2);

end


%-------------------------------------------------------------------------------
function [t, hTE, dtdx, tcp, tdp] = TE_info(X)
% returns trailing-edge information for an airfoil with node coords X
% INPUT
%   X : node coordinates, ordered clockwise (2xN)
% OUTPUT
%   t    : bisector vector = average of upper/lower tangents, normalized
%   hTE  : trailing edge gap, measured as a cross-section
%   dtdx : thickness slope = d(thickness)/d(wake x)
%   tcp  : |t cross p|, used for setting TE source panel strength
%   tdp  : t dot p, used for setting TE vortex panel strength
% DETAILS
%   p refers to the unit vector along the TE panel (from lower to upper)

  t1 = X(:,  1)-X(:,    2); t1 = t1/norm(t1); % lower tangent vector
  t2 = X(:,end)-X(:,end-1); t2 = t2/norm(t2); % upper tangent vector
  t = 0.5*(t1+t2); t = t/norm(t); % average tangent; gap bisector
  s = X(:,end)-X(:,1); % lower to upper connector vector
  hTE = -s(1)*t(2) + s(2)*t(1); % TE gap
  dtdx = t1(1)*t2(2) - t2(1)*t1(2); % sin(theta between t1,t2) approx dt/dx
  p = s/norm(s); % unit vector along TE panel
  tcp = abs(t(1)*p(2)-t(2)*p(1)); tdp = dot(t,p);

end


%-------------------------------------------------------------------------------
function [a,b] = panel_linvortex_velocity(Xj, xi, vdir, onmid)
% calculates the velocity coefficients for a linear vortex panel
% INPUTS
%   Xj    : X(:,[1,2]) = panel endpoint coordinates
%   xi    : control point coordinates (2x1)
%   vdir  : direction of dot product
%   onmid : true means xi is on the panel midpoint
% OUTPUTS
%   a,b   : velocity influence coefficients of the panel
% DETAILS
%   The velocity due to the panel is then a*g1 + b*g2
%   where g1 and g2 are the vortex strengths at the panel endpoints
%   If vdir is [], a,b are 2x1 vectors with velocity components
%   Otherwise, a,b are dotted with vdir
  
  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point
  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle

  % velocity in panel-aligned coord system
  if (onmid)
    ug1 = 1/2 - 1/4;   ug2 = 1/4;
    wg1 = -1/(2*pi);   wg2 = 1/(2*pi);
  else
    temp1 = (theta2-theta1)/(2*pi);
    temp2 = (2*z*log(r1/r2) - 2*x*(theta2-theta1))/(4*pi*d);
    ug1 =  temp1 + temp2;
    ug2 =        - temp2;
    temp1 = log(r2/r1)/(2*pi);
    temp2 = (x*log(r1/r2) - d + z*(theta2-theta1))/(2*pi*d);
    wg1 =  temp1 + temp2;
    wg2 =        - temp2;  
  end

  % velocity influence in original coord system
  a = [ug1*t(1)+wg1*n(1); ug1*t(2)+wg1*n(2)]; % point 1
  b = [ug2*t(1)+wg2*n(1); ug2*t(2)+wg2*n(2)]; % point 2
  if (~isempty(vdir)), a = a'*vdir; b = b'*vdir; end

end


%-------------------------------------------------------------------------------
function [a, b] = panel_linvortex_stream(Xj, xi)
% calculates the streamfunction coefficients for a linear vortex panel
% INPUTS
%   Xj  : X(:,[1,2]) = panel endpoint coordinates
%   xi  : control point coordinates (2x1)
% OUTPUTS
%   a,b : streamfunction influence coefficients
% DETAILS
%   The streamfunction due to the panel is then a*g1 + b*g2
%   where g1 and g2 are the vortex strengths at the panel endpoints

  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point
  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle

  % check for r1, r2 zero
  ep = 1e-10;
  if (r1 < ep), logr1 = 0; else, logr1 = log(r1); end
  if (r2 < ep), logr2 = 0; else, logr2 = log(r2); end
  
  % streamfunction components
  P1 = (0.5/pi)*(z*(theta2-theta1) - d + x*logr1 - (x-d)*logr2);
  P2 = x*P1 + (0.5/pi)*(0.5*r2^2*logr2 - 0.5*r1^2*logr1 - r2^2/4 + r1^2/4);
  
  % influence coefficients
  a = P1-P2/d;
  b =    P2/d;
  
end


%-------------------------------------------------------------------------------
function [a] = panel_constsource_velocity(Xj, xi, vdir)
% calculates the velocity coefficient for a constant source panel
% INPUTS
%   Xj    : X(:,[1,2]) = panel endpoint coordinates
%   xi    : control point coordinates (2x1)
%   vdir  : direction of dot product
% OUTPUTS
%   a     : velocity influence coefficient of the panel
% DETAILS
%   The velocity due to the panel is then a*s
%   where s is the panel source strength
%   If vdir is [], a,b are 2x1 vectors with velocity components
%   Otherwise, a,b are dotted with vdir

  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point
  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle
  
  ep = 1e-9;
  if (r1 < ep), logr1 = 0; theta1=pi; theta2=pi; else, logr1 = log(r1); end
  if (r2 < ep), logr2 = 0; theta1=0; theta2=0; else, logr2 = log(r2); end
  

  % velocity in panel-aligned coord system
  u = (0.5/pi)*(logr1 - logr2);
  w = (0.5/pi)*(theta2-theta1);
  
  % velocity in original coord system dotted with given vector
  a = [u*t(1)+w*n(1); u*t(2)+w*n(2)];
  if (~isempty(vdir)), a = a'*vdir; end

end


%-------------------------------------------------------------------------------
function [a] = panel_constsource_stream(Xj, xi)
% calculates the streamfunction coefficient for a constant source panel
% INPUTS
%   Xj    : X(:,[1,2]) = panel endpoint coordinates
%   xi    : control point coordinates (2x1)
% OUTPUTS
%   a     : streamfunction influence coefficient of the panel
% DETAILS
%   The streamfunction due to the panel is then a*s
%   where s is the panel source strength
  
  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point

  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle
    
  % streamfunction
  ep = 1e-9;
  if (r1 < ep), logr1 = 0; theta1=pi; theta2=pi; else, logr1 = log(r1); end
  if (r2 < ep), logr2 = 0; theta1=0; theta2=0; else, logr2 = log(r2); end
  P = (x*(theta1-theta2) + d*theta2 + z*logr1 - z*logr2)/(2*pi);
  
  dP = d; % delta psi
  if ((theta1+theta2) > pi) 
    P = P - 0.25*dP; 
  else
    P = P + 0.75*dP;
  end

  % influence coefficient
  a = P;
  
end


%-------------------------------------------------------------------------------
function [a, b] = panel_linsource_velocity(Xj, xi, vdir)
% calculates the velocity coefficients for a linear source panel
% INPUTS
%   Xj    : X(:,[1,2]) = panel endpoint coordinates
%   xi    : control point coordinates (2x1)
%   vdir  : direction of dot product
% OUTPUTS
%   a,b   : velocity influence coefficients of the panel
% DETAILS
%   The velocity due to the panel is then a*s1 + b*s2
%   where s1 and s2 are the source strengths at the panel endpoints
%   If vdir is [], a,b are 2x1 vectors with velocity components
%   Otherwise, a,b are dotted with vdir

  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point
  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle

  % velocity in panel-aligned coord system
  temp1 = log(r1/r2)/(2*pi);
  temp2 = (x*log(r1/r2) - d + z*(theta2-theta1))/(2*pi*d);
  ug1 =  temp1 - temp2;
  ug2 =          temp2;
  temp1 = (theta2-theta1)/(2*pi);
  temp2 = (-z*log(r1/r2) + x*(theta2-theta1))/(2*pi*d);
  wg1 =  temp1 - temp2;
  wg2 =          temp2;
  
  % velocity influence in original coord system
  a = [ug1*t(1)+wg1*n(1); ug1*t(2)+wg1*n(2)]; % point 1
  b = [ug2*t(1)+wg2*n(1); ug2*t(2)+wg2*n(2)]; % point 2
  if (~isempty(vdir)), a = a'*vdir; b = b'*vdir; end
  
end


%-------------------------------------------------------------------------------
function [a, b] = panel_linsource_stream(Xj, xi)
% calculates the streamfunction coefficients for a linear source panel
% INPUTS
%   Xj  : X(:,[1,2]) = panel endpoint coordinates
%   xi  : control point coordinates (2x1)
% OUTPUTS
%   a,b : streamfunction influence coefficients
% DETAILS
%   The streamfunction due to the panel is then a*s1 + b*s2
%   where s1 and s2 are the source strengths at the panel endpoints

  % panel coordinates
  xj1 = Xj(1,1);  zj1 = Xj(2,1);
  xj2 = Xj(1,2);  zj2 = Xj(2,2);

  % panel-aligned tangent and normal vectors
  t = [xj2-xj1; zj2-zj1]; t = t/norm(t);
  n = [-t(2); t(1)];

  % control point relative to (xj1,zj1)
  xz = [(xi(1)-xj1), (xi(2)-zj1)];
  x = xz*t;  % in panel-aligned coord system
  z = xz*n;  % in panel-aligned coord system

  % distances and angles
  d = norm([xj2-xj1, zj2-zj1]); % panel length
  r1 = norm([x,z]);             % left edge to control point
  r2 = norm([x-d,z]);           % right edge to control point
  theta1 = atan2(z,x);          % left angle
  theta2 = atan2(z,x-d);        % right angle

  % make branch cut at theta = 0
  if (theta1<0), theta1 = theta1 + 2*pi; end
  if (theta2<0), theta2 = theta2 + 2*pi; end
  
  % check for r1, r2 zero
  ep = 1e-9;
  if (r1 < ep), logr1 = 0; theta1=pi; theta2=pi; else, logr1 = log(r1); end
  if (r2 < ep), logr2 = 0; theta1=0; theta2=0; else, logr2 = log(r2); end

  % streamfunction components
  P1 = (0.5/pi)*(x*(theta1-theta2) + theta2*d + z*logr1 - z*logr2);
  P2 = x*P1 + (0.5/pi)*(0.5*r2^2*theta2 - 0.5*r1^2*theta1 - 0.5*z*d);
  
  % influence coefficients
  a = P1-P2/d;
  b =    P2/d;
  
end


%% ============ GEOMETRY ==============


%-------------------------------------------------------------------------------
function mgeom_flap(M, xzhinge, eta)
% deploys a flap at hinge location xzhinge, with flap angle eta 
% INPUTS
%   M       : mfoil class containing an airfoil
%   xzhinge : flap hinge location (x,z)
%   eta     : flap angle, positive = down, degrees
% OUTPUTS
%   M.foil.x : modified airfoil coordinates

  if (~iscolumn(xzhinge)), xzhinge = xzhinge'; end
  X = M.geom.xpoint; N = size(X,2); % airfoil points
  xh = xzhinge(1);  % x hinge location

  % identify points on flap
  If = find(X(1,:)>xh);

  % rotate flap points
  R = [cosd(eta), sind(eta); -sind(eta), cosd(eta)];
  X(:,If) = xzhinge + R*(X(:,If)-xzhinge);
  
  % remove flap points to left of hinge
  I = If(X(1,If)<xh); I = setdiff(1:N,I);
  
  % re-assemble the airfoil; note, chord length is *not* redefined
  M.geom.xpoint = X(:,I); M.geom.npoint = size(M.geom.xpoint,2); 
    
  % repanel
  if (M.foil.N > 0), make_panels(M, M.foil.N-1); end
  
  % clear the solution
  clear_solution(M);

end


%-------------------------------------------------------------------------------
function mgeom_addcamber(M, xzcamb)
% adds camber to airfoil from given coordinates
% INPUTS
%   M       : mfoil class containing an airfoil
%   xzcamb  : (x,z) points on camberline increment, 2 x Nc
% OUTPUTS
%   M.foil.x : modified airfoil coordinates

  if (size(xzcamb,1) > size(xzcamb,2)), xzcamb = xzcamb'; end
  X = M.geom.xpoint; % airfoil points

  % interpolate camber delta, add to X
  dz = interp1(xzcamb(1,:), xzcamb(2,:), X(1,:), 'spline');
  X(2,:) = X(2,:) + dz;
    
  % store back in M.geom
  M.geom.xpoint = X; M.geom.npoint = size(M.geom.xpoint,2); 
    
  % repanel
  if (M.foil.N > 0), make_panels(M, M.foil.N-1); end
  
  % clear the solution
  clear_solution(M);

end


%-------------------------------------------------------------------------------
function mgeom_derotate(M)
% derotates airfoil about leading edge to make chordline horizontal
% INPUTS
%   M       : mfoil class containing an airfoil
% OUTPUTS
%   M.foil.x : modified airfoil coordinates

  X = M.geom.xpoint; N = size(X,2); % airfoil points
  
  [~,I] = min(X(1,:)); xLE = X(:,I(1)); % LE point
  xTE = 0.5*(X(:,1) + X(:,N)); % TE "point"
  
  theta = atan2(xTE(2)-xLE(2), xTE(1)-xLE(1)); % rotation angle
  R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
  X = xLE + R*(X-xLE); % rotation
    
  % store back in M.geom
  M.geom.xpoint = X; M.geom.npoint = size(M.geom.xpoint,2); 
    
  % repanel
  if (M.foil.N > 0), make_panels(M, M.foil.N-1); end
  
  % clear the solution
  clear_solution(M);

end


%-------------------------------------------------------------------------------
function x = space_geom(dx0, L, Np)
% spaces Np points geometrically from [0,L], with dx0 as first interval
% INPUTS
%   dx0 : first interval length
%   L   : total domain length
%   Np  : number of points, including endpoints at 0,L
% OUTPUTS
%   x   : point locations (1xN)
    
  assert(Np>1, 'Need at least two points for spacing.');
  N = Np - 1; % number of intervals
  % L = dx0 * (1 + r + r^2 + ... r^{N-1}) = dx0*(r^N-1)/(r-1)
  % let d = L/dx0, and for a guess, consider r = 1 + s
  % The equation to solve becomes d*s  = (1+s)^N - 1
  % Initial guess: (1+s)^N ~ 1 + N*s + N*(N-1)*s^2/2 + N*(N-1)*(N-2)*s^3/3
  d = L/dx0; a = N*(N-1.)*(N-2.)/6.; b = N*(N-1.)/2.; c = N-d;
  disc = max(b*b-4.*a*c, 0.); r = 1 + (-b+sqrt(disc))/(2*a);
  for k = 1:10
    R = r^N -1-d*(r-1); R_r = N*r^(N-1)-d; dr = -R/R_r;
    if (abs(dr)<1e-6), break; end; r = r - R/R_r;
  end
  x = [0,cumsum(dx0*r.^(0:N-1))];
  
end


%-------------------------------------------------------------------------------
function set_coords(M, X)
% sets geometry from coordinate matrix
% INPUTS
%   M : mfoil class
%   X : matrix whose rows or columns are (x,z) points, CW or CCW
% OUTPUTS
%   M.geom.npoint : number of points
%   M.geom.xpoint : point coordinates (2 x npoint)
%   M.geom.chord  : chord length
% DETAILS
%   Coordinates should start and end at the trailing edge
%   Trailing-edge point must be repeated if sharp
%   Points can be clockwise or counter-clockwise (will detect and make CW)

  if (size(X,1) > size(X,2)), X = X'; end

  % ensure CCW
  A = 0.;
  for i = 1:size(X,2)-1, A = A + (X(1,i+1)-X(1,i))*(X(2,i+1)+X(2,i)); end
  if (A<0), X = fliplr(X); end
  
  % store points in M
  M.geom.npoint = size(X,2);
  M.geom.xpoint = X;  
  M.geom.chord = max(X(1,:)) - min(X(1,:));

end


%-------------------------------------------------------------------------------
function naca_points(M, digits)
% calculates coordinates of a NACA 4-digit airfoil, stores in M.geom
% INPUTS
%   M      : mfoil class
%   digits : character array containing NACA digits
% OUTPUTS
%   M.geom.npoint : number of points
%   M.geom.xpoint : point coordinates (2 x npoint)
%   M.geom.chord  : chord length
% DETAILS
%   Uses analytical camber/thickness formulas

  M.geom.name = sprintf('NACA %s', digits);
  N = 100; te = 1.5; % # points per side and trailing-edge bunching factor
  f = linspace(0,1,N+1); % linearly-spaced points between 0 and 1
  x = 1 - (te+1)*f.*(1-f).^te - (1-f).^(te+1); % bunched points, x, 0 to 1
  % normalized thickness, gap at trailing edge (use -.1036*x^4 for no gap)
  t = .2969*sqrt(x) - .126*x - .3516*x.^2 + .2843*x.^3 - .1015*x.^4;
  tmax = str2double(digits(end-1:end))*0.01; % maximum thickness
  t = t*tmax/.2;
  
  
  if (length(digits)==4)
    % 4-digit series
    m = str2double(digits(1))*0.01;
    p = str2double(digits(2))*0.1;
    c = m/(1-p)^2 * ((1-2.*p)+2.*p*x-x.^2);
    I = find(x<p); c(I) = m/p^2*(2*p*x(I)-x(I).^2);
  elseif (length(digits)==5)
    % 5-digit series
    n = str2double(digits(2)); valid = prod(digits([1,3])=='20') && (n>0) && (n<6);
    assert(valid, '5-digit NACA must begin with 2X0, X in 1-5');
    mv = [.058, .126, .2025, .29, .391]; m = mv(n);
    cv = [361.4, 51.64, 15.957, 6.643, 3.23]; cc = cv(n);
    c = (cc/6)*(x.^3 - 3*m*x.^2 + m^2*(3-m)*x);
    I = find(x>m); c(I) = (cc/6)*m^3*(1-x(I));
  else
    error('Provide 4 or 5 NACA digits');
  end

  zu = c + t; zl = c-t; % upper and lower surfaces
  xs = [fliplr(x), x(2:end)]; % x points
  zs = [fliplr(zl), zu(2:end)];% z points

  % store points in M
  M.geom.npoint = length(xs);
  M.geom.xpoint = [xs; zs];  
  M.geom.chord = max(xs) - min(xs);
  
end


%-------------------------------------------------------------------------------
function [X, S, XS] = spline_curvature(Xin, N, Ufac, TEfac)
% Splines 2D points in Xin and samples using curvature-based spacing 
% INPUT
%   Xin   : points to spline
%   N     : number of points = one more than the number of panels
%   Ufac  : uniformity factor (1 = normal; > 1 means more uniform distribution)
%   TEfac : trailing-edge resolution factor (1 = normal; > 1 = high; < 1 = low)
% OUTPUT
%   X  : new points (2xN)
%   S  : spline s values (N)
%   XS : spline tangents (2xN)
  
  % min/max of given points (x-coordinate)
  xmin = min(Xin(1,:)); xmax = max(Xin(1,:));
  
  % spline given points
  PP = spline2d(Xin);

  % curvature-based spacing on geom
  nfine = 501;
  s = linspace(0,PP.breaks(end),nfine);
  xyfine = splineval(PP, s);
  PPfine = spline2d(xyfine);
  s = PPfine.breaks;
  
  sk = zeros(1,nfine);
  [xq, wq] = quadseg;
  for i = 1:nfine-1
    ds = s(i+1)-s(i);
    st = xq*ds;
    px = PPfine.xcoefs(i,:);
    xss = 6.0*px(1)*st + 2.0*px(2);
    py = PPfine.ycoefs(i,:);
    yss = 6.0*py(1)*st + 2.0*py(2);
    skint = 0.01*Ufac+0.5*dot(wq, sqrt(xss.*xss + yss.*yss))*ds;
  
    % force TE resolution
    xx = (0.5*(xyfine(1,i)+xyfine(1,i+1))-xmin)/(xmax-xmin); % close to 1 means at TE
    skint = skint + TEfac*0.5*exp(-100*(1.0-xx));

    % increment sk
    sk(i+1) = sk(i) + skint;
  end
  
  % offset by fraction of average to avoid problems with zero curvature
  sk = sk + 2.0*sum(sk)/nfine;
 
  % arclength values at points
  skl = linspace(min(sk), max(sk), N);
  s = interp1(sk, s, skl, 'spline');
 
  % new points
  X  = splineval(PPfine, s); S = s; XS = splinetan(PPfine, s);

end


%-------------------------------------------------------------------------------
function [PP] = spline2d(X)
% splines 2d points
% INPUT
%   X : points to spline (2xN)
% OUTPUT
%   PP : two-dimensional spline structure 
  
  N = length(X); S = zeros(1,N); Snew = zeros(1,N);

  % estimate the arclength and spline x, y separately
  for i=2:N, S(i) = S(i-1) + norm(X(:,i)-X(:,i-1)); end
  PPX = spline(S,X(1,:)); PPY = spline(S,X(2,:));
  
  % re-integrate to true arclength via several passes
  [xq, wq] = quadseg;
  for ipass = 1:10
    serr = 0;
    Snew(1) = S(1);
    for i = 1:(N-1)
      ds = S(i+1)-S(i);
      st = xq*ds;
      px = PPX.coefs(i,:);
      xs = 3.0*px(1)*st.^2 + 2.0*px(2)*st + px(3);
      py = PPY.coefs(i,:);
      ys = 3.0*py(1)*st.^2 + 2.0*py(2)*st + py(3);
      sint = dot(wq, sqrt(xs.*xs + ys.*ys))*ds;
      serr = max(serr, abs(sint-ds));
      Snew(i+1) = Snew(i) + sint;
    end
    S = Snew;
    PPX = spline(S,X(1,:));
    PPY = spline(S,X(2,:));
  end
  
  PP = splinetwo2one(PPX, PPY);

end
  
  
%-------------------------------------------------------------------------------
function [XY] = splineval(PP, S)
% evaluates 2d spline at given S values
% INPUT
%   PP : two-dimensional spline structure 
%   S  : arclength values at which to evaluate the spline
% OUTPUT
%   XY : coordinates on spline at the requested s values (2xN)
  
  [PPX, PPY] = splineone2two(PP);
  XY = [ppval(PPX, S); ppval(PPY, S)];
  
end


%-------------------------------------------------------------------------------
function [XYS] = splinetan(PP, S)
% evaluates 2d spline tangent (not normalized) at given S values
% INPUT
%   PP  : two-dimensional spline structure 
%   S   : arclength values at which to evaluate the spline tangent
% OUTPUT
%   XYS : dX/dS and dY/dS values at each point (2xN)
  
  [PPX, PPY] = splineone2two(PP);
  C = [diag([3,2,1]); zeros(1,3)];
  PPX.coefs = PPX.coefs*C; PPX.order = 3;
  PPY.coefs = PPY.coefs*C; PPY.order = 3;
  XYS = [ppval(PPX,S); ppval(PPY,S)];  
  
end
  

%-------------------------------------------------------------------------------
function [PP] = splinetwo2one(PPX, PPY)
% combines separate x,y splines into one 2d spline
% INPUT
%   PPX, PPY : one-dimensional spline structures
% OUTPUT
%   PP : two-dimensional spline structure
  
  c = ((PPX.pieces==PPY.pieces) & (PPX.order==PPY.order) & (PPX.dim==PPY.dim));
  assert(c, 'Splines not compatible in 1d->2d merge.')
  PP.breaks  = PPX.breaks;
  PP.xcoefs  = PPX.coefs;
  PP.ycoefs  = PPY.coefs;
  PP.pieces  = PPX.pieces;
  PP.order   = PPX.order;
  PP.dim     = PPX.dim;
  
end
  

%-------------------------------------------------------------------------------
function [PPX, PPY] = splineone2two(PP)
% splits a 2d spline into two 1d splines
% INPUT
%   PP : two-dimensional spline structure
% OUTPUT
%   PPX, PPY : one-dimensional spline structures

  PPX = mkpp(PP.breaks, PP.xcoefs);
  PPY = mkpp(PP.breaks, PP.ycoefs);
  
end


%-------------------------------------------------------------------------------
function [x, w] = quadseg
% Returns quadrature points and weights for a [0,1] line segment
% INPUT
% OUTPUT
%   x : quadrature point coordinates (1d)
%   w : quadrature weights
    
  x = [ 0.046910077030668, 0.230765344947158, 0.500000000000000, ...
        0.769234655052842, 0.953089922969332];
  w = [ 0.118463442528095, 0.239314335249683, 0.284444444444444, ...
        0.239314335249683, 0.118463442528095];
end 


%% ============ VISCOUS FUNCTIONS ==============


%-------------------------------------------------------------------------------
function calc_ue_m(M)
% calculates sensitivity matrix of ue w.r.t. transpiration BC mass sources
% INPUT
%   M : mfoil class with wake already built
% OUTPUT
%   M.vsol.sigma_m : d(source)/d(mass) matrix, for computing source strengths
%   M.vsol.ue_m    : d(ue)/d(mass) matrix, for computing tangential velocity
% DETAILS
%   "mass" flow refers to area flow (we exclude density)
%   sigma_m and ue_m return values at each node (airfoil and wake)
%   airfoil panel sources are constant strength
%   wake panel sources are two-piece linear
  
  assert(~isempty(M.isol.gam), 'No inviscid solution');
  N = M.foil.N; Nw = M.wake.N;  % number of points on the airfoil/wake
  assert(Nw>0, 'No wake');
  
  % Cgam = d(wake uei)/d(gamma)   [Nw x N]   (not sparse)
  Cgam = zeros(Nw,N); 
  for i = 1:Nw
    [~, v_G] = inviscid_velocity(M.foil.x, M.isol.gam, 0, 0, M.wake.x(:,i));
    Cgam(i,:) = v_G(1,:)*M.wake.t(1,i) + v_G(2,:)*M.wake.t(2,i);
  end
  
  % B = d(airfoil surf streamfunction)/d(source)  [(N+1) x (N+Nw-2)]  (not sparse)
  B = zeros(N+1,N+Nw-2);  % note, N+Nw-2 = # of panels
  for i = 1:N  % loop over points on the airfoil
    xi = M.foil.x(:,i); % coord of point i
    for j = 1:N-1 % loop over airfoil panels
      B(i,j) = panel_constsource_stream(M.foil.x(:,[j,j+1]), xi);
    end
    for j = 1:Nw-1 % loop over wake panels
      Xj = M.wake.x(:,[j,j+1]); % panel endpoint coordinates
      Xm = 0.5*(Xj(:,1) + Xj(:,2)); % panel midpoint
      Xj = [Xj(:,1), Xm, Xj(:,2)]; % left, mid, right coords on panel
      if (j==(Nw-1)), Xj(:,3) = 2*Xj(:,3) - Xj(:,2); end % ghost extension at last point
      [a,b] = panel_linsource_stream(Xj(:,[1,2]), xi); % left half panel
      if (j > 1)
        B(i,N-1+j) = B(i,N-1+j) + 0.5*a + b;
        B(i,N-1+j-1) = B(i,N-1+j-1) + 0.5*a;
      else
        B(i,N-1+j) = B(i,N-1+j) + b;
      end
      [a,b] = panel_linsource_stream(Xj(:,[2,3]), xi); % right half panel
      B(i,N-1+j) = B(i,N-1+j) + a + 0.5*b;
      if (j<Nw-1)
        B(i,N-1+j+1) = B(i,N-1+j+1) + 0.5*b;
      else
        B(i,N-1+j) = B(i,N-1+j) + 0.5*b;
      end
    end
  end
  
  % Bp = - inv(AIC) * B   [N x (N+Nw-2)]  (not sparse)
  % Note, Bp is d(airfoil gamma)/d(source)
  Bp = -M.isol.AIC\B;  % this has N+1 rows, but the last one is zero
  Bp = Bp(1:end-1,:);  % trim the last row
  
  % Csig = d(wake uei)/d(source) [Nw x (N+Nw-2)]  (not sparse)
  Csig = zeros(Nw, N+Nw-2);
  for i = 1:Nw
    xi = M.wake.x(:,i); ti = M.wake.t(:,i); % point, tangent on wake
    
    % first/last airfoil panel effects on i=1 wake point handled separately
    jstart = 1 + (i==1); jend = N-1 - (i==1);
    for j = jstart:jend % constant sources on airfoil panels
      Csig(i,j) = panel_constsource_velocity(M.foil.x(:,[j,j+1]), xi, ti);
    end
    
    % piecewise linear sources across wake panel halves (else singular)
    for j = 1:Nw  % loop over wake points
      I = [max(j-1,1), j, min(j+1,Nw)]; % left, self, right
      Xj = M.wake.x(:,I); % point coordinates
      Xj(:,1) = 0.5*(Xj(:,1) + Xj(:,2)); % left midpoint
      Xj(:,3) = 0.5*(Xj(:,2) + Xj(:,3)); % right midpoint
      if (j==Nw), Xj(:,3) = 2*Xj(:,2) - Xj(:,1); end % ghost extension at last point
      d1 = norm(Xj(:,2)-Xj(:,1)); % left half-panel length
      d2 = norm(Xj(:,3)-Xj(:,2)); % right half-panel length
      if (i==j)
        if (j==1) % first point: special TE system (three panels meet)
          dl = norm(M.foil.x(:,2)-M.foil.x(:,  1)); % lower surface panel length
          du = norm(M.foil.x(:,N)-M.foil.x(:,N-1)); % upper surface panel length
          Csig(i,  1) = Csig(i,  1) + (0.5/pi)*(log(dl/d2) + 1); % lower panel effect
          Csig(i,N-1) = Csig(i,N-1) + (0.5/pi)*(log(du/d2) + 1); % upper panel effect
          Csig(i,N-1+1) = Csig(i,N-1+1) - 0.5/pi; % self effect
        elseif (j==Nw) % last point: no self effect of last pan (ghost extension)
          Csig(i,N-1+j-1) = Csig(i,N-1+j-1) + 0; % hence the 0
        else % all other points
          aa = (0.25/pi)*log(d1/d2);
          Csig(i,N-1+j-1) = Csig(i,N-1+j-1) + aa + 0.5/pi;
          Csig(i,N-1+j  ) = Csig(i,N-1+j  ) + aa - 0.5/pi;
        end
      else
        if (j==1) % first point only has a half panel on the right
          [a,b] = panel_linsource_velocity(Xj(:,[2,3]), xi, ti);
          Csig(i,N-1+1) = Csig(i,N-1+1) + b; % right half panel effect
          Csig(i,  1  ) = Csig(i,  1  ) + a; % lower airfoil panel effect
          Csig(i,N-1  ) = Csig(i,N-1  ) + a; % upper airfoil panel effect
        elseif (j==Nw) % last point has a constant source ghost extension
          a = panel_constsource_velocity(Xj(:,[1,3]), xi, ti);
          Csig(i,N+Nw-2) = Csig(i,N+Nw-2) + a; % full const source panel effect
        else % all other points have a half panel on left and right
          [a1,b1] = panel_linsource_velocity(Xj(:,[1,2]), xi, ti); % left half-panel ue contrib
          [a2,b2] = panel_linsource_velocity(Xj(:,[2,3]), xi, ti); % right half-panel ue contrib
          Csig(i,N-1+j-1) = Csig(i,N-1+j-1) + a1 + 0.5*b1;
          Csig(i,N-1+j  ) = Csig(i,N-1+j  ) + 0.5*a2 + b2;
        end
      end
    end
  end
  
  % compute ue_sigma = d(unsigned ue)/d(source) [(N+Nw) x (N+Nw-2)] (not sparse)
  % Df = +/- Bp = d(foil uei)/d(source)  [N x (N+Nw-2)]  (not sparse)
  % Dw = (Cgam*Bp + Csig) = d(wake uei)/d(source)  [Nw x (N+Nw-2)]  (not sparse)
  Dw = Cgam*Bp + Csig;
  Dw(1,:) = Bp(end,:); % ensure first wake point has same ue as TE
  M.vsol.ue_sigma = [Bp; Dw]; % store combined matrix
  
  % build ue_m from ue_sigma, using sgnue
  rebuild_ue_m(M);

end


%-------------------------------------------------------------------------------
function rebuild_ue_m(M)
% rebuilds ue_m matrix after stagnation panel change (new sgnue)
% INPUT
%   M : mfoil class with calc_ue_m already called once
% OUTPUT
%   M.vsol.sigma_m : d(source)/d(mass) matrix, for computing source strengths
%   M.vsol.ue_m    : d(ue)/d(mass) matrix, for computing tangential velocity
% DETAILS
%   "mass" flow refers to area flow (we exclude density)
%   sigma_m and ue_m return values at each node (airfoil and wake)
%   airfoil panel sources are constant strength
%   wake panel sources are two-piece linear
  
  assert(~isempty(M.vsol.ue_sigma), 'Need ue_sigma to build ue_m');
  
  % Dp = d(source)/d(mass)  [(N+Nw-2) x (N+Nw)]  (sparse)
  N = M.foil.N; Nw = M.wake.N;  % number of points on the airfoil/wake
  Dp = spalloc(N+Nw-2,N+Nw, 2*(N-1+Nw-1));
  for i = 1:N-1
    ds = M.foil.s(i+1)-M.foil.s(i);
    % Note, at stagnation: ue = K*s, dstar = const, m = K*s*dstar
    % sigma = dm/ds = K*dstar = m/s (separate for each side, +/-)
    Dp(i,[i,i+1]) = M.isol.sgnue([i,i+1]).*[-1,1]/ds;
  end
  for i = 1:Nw-1
    ds = M.wake.s(i+1)-M.wake.s(i);
    Dp(N-1+i,[N+i,N+i+1]) = [-1,1]/ds;
  end
  M.vsol.sigma_m = Dp;
  
  % sign of ue at all points (wake too)
  sgue = [M.isol.sgnue, ones(1,Nw)];
  
  % ue_m = ue_sigma * sigma_m [(N+Nw) x (N+Nw)] (not sparse)
  M.vsol.ue_m = spdiags(sgue',0,N+Nw,N+Nw)*M.vsol.ue_sigma*M.vsol.sigma_m;
  
end


%-------------------------------------------------------------------------------
function init_thermo(M)
% initializes thermodynamics variables in param structure
% INPUT
%   M  : mfoil class with oper structure set
% OUTPUT
%   M.param fields filled in based on M.oper
%   Gets ready for compressibilty corrections if M.oper.Ma > 0

  g = M.param.gam; gmi = g-1;
  rhoinf = M.oper.rho; % freestream density
  Vinf = M.oper.Vinf; M.param.Vinf = Vinf; % freestream speed
  M.param.muinf = rhoinf*Vinf*M.geom.chord/M.oper.Re; % freestream dyn viscosity 
  Minf = M.oper.Ma; M.param.Minf = Minf; % freestream Mach
  if (Minf > 0)
    M.param.KTb = sqrt(1-Minf^2); % Karman-Tsien beta
    M.param.KTl = Minf^2/(1+M.param.KTb)^2; % Karman-Tsien lambda
    M.param.H0 = (1+0.5*gmi*Minf^2)*Vinf^2/(gmi*Minf^2); % stagnation enthalpy
    Tr = 1-0.5*Vinf^2/M.param.H0; % freestream/stagnation temperature ratio
    finf = Tr^1.5*(1+M.param.Tsrat)/(Tr + M.param.Tsrat); % Sutherland's ratio
    M.param.cps = 2/(g*Minf^2)*(((1+0.5*gmi*Minf^2)/(1+0.5*gmi))^(g/gmi) - 1);
  else
    finf = 1; % incompressible case
  end
  M.param.mu0 = M.param.muinf/finf;  % stag visc (Sutherland ref temp is stag)
  M.param.rho0 = rhoinf*(1+0.5*gmi*Minf^2)^(1/gmi); % stag density

end


%-------------------------------------------------------------------------------
function identify_surfaces(M)
% identifies lower/upper/wake surfaces
% INPUT
%   M  : mfoil class with stagnation point found
% OUTPUT
%   M.vsol.Is : cell array of node indices for lower(1), upper(2), wake(3)

  M.vsol.Is{1} = M.isol.Istag(1):-1:1;
  M.vsol.Is{2} = M.isol.Istag(2):M.foil.N;
  M.vsol.Is{3} = (M.foil.N+1):(M.foil.N+M.wake.N);

end


%-------------------------------------------------------------------------------
function set_wake_gap(M)
% sets height (delta*) of dead air in wake
% INPUT
%   M  : mfoil class with wake built and stagnation point found
% OUTPUT
%   M.vsol.wgap : wake gap at each wake point
% DETAILS
%   Uses cubic function to extrapolate the TE gap into the wake
%   See Drela, IBL for Blunt Trailing Edges, 1989, 89-2166-CP
  
  [~, hTE, dtdx,~,~] = TE_info(M.foil.x);
  flen = 2.5; % length-scale factor
  dtdx = min(max(dtdx,-3./flen), 3./flen); % clip TE thickness slope
  Lw = flen*hTE;
  wgap = zeros(1,M.wake.N);
  for i = 1:M.wake.N
    xib = (M.isol.xi(M.foil.N+i) - M.isol.xi(M.foil.N))/Lw;
    if (xib <= 1), wgap(i) = hTE*(1+(2+flen*dtdx)*xib)*(1-xib)^2; end
  end
  M.vsol.wgap = wgap;
  
end


%-------------------------------------------------------------------------------
function stagpoint_move(M)
% moves the LE stagnation point on the airfoil using the global solution ue
% INPUT
%   M  : mfoil class with a valid solution in M.glob.U
% OUTPUT
%   New sstag, sstag_ue, xi in M.isol
%   Possibly new stagnation panel, Istag, and hence new surfaces and matrices
  
  N = M.foil.N;  % number of points on the airfoil
  I = M.isol.Istag; % current adjacent node indices
  ue = M.glob.U(4,:)'; % edge velocity
  sstag0 = M.isol.sstag; % original stag point location
  
  newpanel = true; % are we moving to a new panel?
  if (ue(I(2)) < 0)
    % move stagnation point up (larger s, new panel)
    vprint(M.param,2, '  Moving stagnation point up\n');
    J = find(ue(I(2):end) > 0); I2 = J(1)+I(2)-1;
    for j = I(2):(I2-1), ue(j) = -ue(j); end
    I = [I2-1, I2]; % new panel
  elseif (ue(I(1)) < 0)
    % move stagnation point down (smaller s, new panel)
    vprint(M.param,2, '  Moving stagnation point down\n');
    J = find(ue(I(1):-1:1) > 0); I1 = I(1)-J(1)+1;
    for j = (I1+1):I(1), ue(j) = -ue(j); end
    I = [I1, I1+1]; % new panel
  else, newpanel = false; % staying on the current panel
  end
  
  % move point along panel
  ues = ue(I); S = M.foil.s(I);
  assert((ues(1) > 0) && (ues(2) > 0), 'stagpoint_move: velocity error');
  den = ues(1) + ues(2); w1 = ues(2)/den; w2 = ues(1)/den;
  M.isol.sstag = w1*S(1) + w2*S(2);  % s location
  M.isol.xstag = M.foil.x(:,I)*[w1;w2];  % x location
  M.isol.sstag_ue = [ues(2), -ues(1)]*(S(2)-S(1))/(den*den);
  vprint(M.param,2, '  Moving stagnation point: s=%.15e -> s=%.15e\n', sstag0, M.isol.sstag);
  
  % set new xi coordinates for every point
  M.isol.xi = [abs(M.foil.s-M.isol.sstag), M.wake.s-M.isol.sstag];
  
  % matrices need to be recalculated if on a new panel
  if (newpanel)
    vprint(M.param,2, '  New stagnation panel = %d %d\n', I(1), I(2));
    M.isol.Istag = I; % new panel indices
    sgnue = ones(1,N); sgnue(1:I(1)) = -1;
    M.isol.sgnue = sgnue; % new upper/lower surface signs
    identify_surfaces(M); % re-identify surfaces
    M.glob.U(4,:) = ue;  % sign of ue changed on some points near stag    
    rebuild_ue_m(M);
  end
    
end

%-------------------------------------------------------------------------------
function solve_viscous(M)
% solves the viscous system (BL + outer flow concurrently)
% INPUT
%   M  : mfoil class with an airfoil
% OUTPUT
%   M.glob.U : global solution
%   M.post   : post-processed quantities

  solve_inviscid(M);
  M.oper.viscous = true;
  init_thermo(M); % thermodynamics
  build_wake(M);
  stagpoint_find(M); % from the inviscid solution
  identify_surfaces(M);
  set_wake_gap(M);  % blunt TE dead air extent in wake
  calc_ue_m(M);
  init_boundary_layer(M); % initialize boundary layer from ue
  stagpoint_move(M); % move stag point, using viscous solution  
  solve_coupled(M); % solve coupled system
  calc_force(M);
  get_distributions(M);
  
end


%-------------------------------------------------------------------------------
function solve_coupled(M)
% Solves the coupled inviscid and viscous system
% INPUT
%   M  : mfoil class with an inviscid solution
% OUTPUT
%   M.glob.U : global coupled solution
% DETAILS
%   Inviscid solution should exist, and BL variables should be initialized
%   The global variables are [th, ds, sa, ue] at every node
%   th = momentum thickness; ds = displacement thickness
%   sa = amplification factor or sqrt(ctau); ue = edge velocity
%   Nsys = N + Nw = total number of unknowns
%   ue is treated as a separate variable for improved solver robustness
%   The alternative is to eliminate ue, ds and use mass flow (not done here):
%     Starting point: ue = uinv + D*m -> ue_m = D
%     Since m = ue*ds, we have ds = m/ue = m/(uinv + D*m)
%     So, ds_m = diag(1/ue) - diag(ds/ue)*D
%     The residual linearization is then: R_m = R_ue*ue_m + R_ds*ds_m
  
  % Newton loop
  nNewton = M.param.niglob; % number of iterations
  M.glob.conv = false;
  vprint(M.param,1, '\n <<< Beginning coupled solver iterations >>> \n');
  for iNewton = 1:nNewton
    
    % set up the global system
    build_glob_sys(M);
    
    % compute forces
    calc_force(M);
    
    % convergence check
    Rnorm = norm(M.glob.R, 2);
    vprint(M.param,1, '\nNewton iteration %d, Rnorm = %.5e\n', iNewton, Rnorm);
    if (Rnorm < M.param.rtol), M.glob.conv = true; break; end
    
    % solve global system
    solve_glob(M);
    
    % update the state
    update_state(M);
        
    % update stagnation point; Newton still OK; had R_x effects in R_U
    stagpoint_move(M);
    
    % update transition
    update_transition(M);
    
  end
  
  if (~M.glob.conv), vprint(M.param,1, '\n** Global Newton NOT CONVERGED **\n'); end
  
end


%-------------------------------------------------------------------------------
function update_state(M)
% updates state, taking into account physical constraints
% INPUT
%   M  : mfoil class with a valid solution (U) and proposed update (dU)
% OUTPUT
%   M.glob.U : updated solution, possibly with a fraction of dU added
% DETAILS
%   U = U + omega * dU; omega = under-relaxation factor
%   Calculates omega to prevent big changes in the state or negative values
    
  if (any(imag(M.glob.U(3,:)))), M.glob.U(3,:), error('imaginary amp in U'); end
  if (any(imag(M.glob.dU(3,:)))), M.glob.dU(3,:), error('imaginary amp in dU'); end

  % max ctau
  It = find(M.vsol.turb); ctmax = max(M.glob.U(3,It));
  
  % starting under-relaxation factor
  omega = 1.0;
  
  % first limit theta and delta*
  for k = 1:2,
    Uk = M.glob.U(k,:); dUk = M.glob.dU(k,:);
    % prevent big decreases in th, ds
    fmin = min(dUk./Uk); % find most negative ratio
    if (fmin < -0.5), om = abs(0.5/fmin); else, om = 1; end
    if (om<omega), omega = om; vprint(M.param,3, '  th/ds decrease: omega = %.5f\n', omega); end
  end
  
  % limit negative amp/ctau
  Uk = M.glob.U(3,:); dUk = M.glob.dU(3,:);
  for i = 1:length(Uk)
    if (~M.vsol.turb(i)) && (Uk(i)<.2), continue; end % do not limit very small amp (too restrictive)
    if (M.vsol.turb(i)) && (Uk(i)<0.1*ctmax), continue; end % do not limit small ctau
    if (Uk(i)==0.) || (dUk(i)==0.), continue; end
    if (Uk(i)+dUk(i) < 0)
      om = 0.8*abs(Uk(i)/dUk(i)); 
      if (om<omega), omega = om; vprint(M.param,3, '  neg sa: omega = %.5f\n', omega); end
    end
  end
  
  % prevent big changes in amp
  I = find(~M.vsol.turb); 
  if (any(imag(Uk(I)))), error('imaginary amplification'); end
  dumax = max(abs(dUk(I)));
  if (dumax > 0), om = abs(2/dumax); else, om = 1; end
  if (om<omega), omega = om; vprint(M.param,3, '  amp: omega = %.5f\n', omega); end
  
  % prevent big changes in ctau
  I = find(M.vsol.turb); 
  dumax = max(abs(dUk(I)));
  if (dumax > 0), om = abs(.05/dumax); else, om = 1; end
  if (om<omega), omega = om; vprint(M.param,3, '  ctau: omega = %.5f\n', omega); end
  
  % prevent large ue changes
  dUk = M.glob.dU(4,:);
  fmax = max(abs(dUk)/M.oper.Vinf);
  if (fmax > 0), om =.2/fmax; else, om = 1; end
  if (om<omega), omega = om; vprint(M.param,3, '  ue: omega = %.5f\n', omega); end
  
  % prevent large alpha changes
  if (abs(M.glob.dalpha) > 2), omega = min(omega, abs(2/M.glob.dalpha)); end

  % take the update
  vprint(M.param,2, '  state update: under-relaxation = %.5f\n', omega);
  M.glob.U = M.glob.U + omega*M.glob.dU;
  M.oper.alpha = M.oper.alpha + omega*M.glob.dalpha;
  
  % fix bad Hk after the update
  for is = 1:3  % loop over surfaces
    if (is==3), Hkmin = 1.00005; else, Hkmin = 1.02; end
    Is = M.vsol.Is{is}; % surface point indices
    param = build_param(M, is); % get parameter structure
    for i = 1:length(Is) % loop over points
      j = Is(i); Uj = M.glob.U(:,j);
      param = station_param(M, param, j);
      [Hk, ~] = get_Hk(Uj, param);
      if (Hk < Hkmin)
        M.glob.U(2,j) = M.glob.U(2,j) + 2*(Hkmin-Hk)*M.glob.U(1,j);
      end
    end
  end
  
  % fix negative ctau after the update
  for ii = 1:length(I)
    i = It(ii); if (M.glob.U(3,i) < 0), M.glob.U(3,i) = 0.1*ctmax; end
  end
  
  % rebuild inviscid solution (gam, wake) if angle of attack changed
  if (abs(omega*M.glob.dalpha) > 1e-10), rebuild_isol(M); end
  
end


%-------------------------------------------------------------------------------
function jacobian_add_Rx(M)
% include effects of R_x into R_U: R_ue += R_x*x_st*st_ue
% INPUT
%   M  : mfoil class with residual Jacobian calculated
% OUTPUT
%   M.glob.R_U : ue linearization updated with R_x
% DETAILS
%   The global residual Jacobian has a column for ue sensitivity
%   ue, the edge velocity, also affects the location of the stagnation point
%   The location of the stagnation point (st) dictates the x value at each node
%   The residual also depends on the x value at each node (R_x)
%   We use the chain rule (formula above) to account for this
  
  Nsys = M.glob.Nsys; % number of dofs
  Iue = 4:4:4*Nsys; % ue indices in U
  x_st = -M.isol.sgnue';  % st = stag point [Nsys x 1]
  x_st = [x_st; -ones(M.wake.N,1)]; % wake same sens as upper surface
  R_st = M.glob.R_x*x_st; % [3*Nsys x 1]
  Ist = M.isol.Istag; st_ue = M.isol.sstag_ue; % stag points, sens
  M.glob.R_U(:,Iue(Ist)) = M.glob.R_U(:,Iue(Ist)) + R_st*st_ue;

end

    
%-------------------------------------------------------------------------------
function solve_glob(M)
% solves global system for the primary variable update dU
% INPUT
%   M  : mfoil class with residual and Jacobian calculated
% OUTPUT
%   M.glob.dU : proposed solution update
% DETAILS
%   Uses the augmented system: fourth residual = ue equation
%   Supports lift-constrained mode, with an extra equation: cl - cltgt = 0
%   Extra variable in cl-constrained mode is angle of attack
%   Solves sparse matrix system for for state/alpha update  
  
  Nsys = M.glob.Nsys; % number of dofs
  docl = M.oper.givencl; % 1 if in cl-constrained mode
  
  % get edge velocity and displacement thickness
  ue = M.glob.U(4,:)'; ds = M.glob.U(2,:)';
  uemax = max(abs(ue)); ue = max(ue,1e-10*uemax); % avoid 0/negative ue
    
  % use augmented system: variables = th, ds, sa, ue
    
  % inviscid edge velocity on the airfoil and wake
  ueinv = get_ueinv(M);
  
  % initialize the global variable Jacobian (TODO: estimate nnz)
  R_V = sparse(4*Nsys+docl,4*Nsys+docl); % +1 for cl-alpha constraint
  
  % state indices in the global system
  Ids = 2:4:4*Nsys; % delta star indices
  Iue = 4:4:4*Nsys; % ue indices
  
  % include effects of R_x into R_U: R_ue += R_x*x_st*st_ue
  jacobian_add_Rx(M);

  % assemble the residual
  R = [M.glob.R; ue - (ueinv + M.vsol.ue_m*(ds.*ue))];
  
  % assemble the Jacobian
  R_V(1:3*Nsys,1:4*Nsys) = M.glob.R_U;
  I = (3*Nsys+1):4*Nsys;
  R_V(I,Iue) = speye(Nsys) - M.vsol.ue_m*diag(ds);
  R_V(I,Ids) = -M.vsol.ue_m*diag(ue);
  
  if (docl)
    % include cl-alpha residual and Jacobian
    [Rcla, Ru_alpha, Rcla_U] = clalpha_residual(M); 
    R = [R; Rcla]; R_V(I,4*Nsys+1) = Ru_alpha; R_V(4*Nsys+1,:) = Rcla_U;
  end
  
  % solve system for dU, dalpha
  dV = -R_V\R;

  % store dU, reshaped, in M
  M.glob.dU = reshape(dV(1:4*Nsys),4,Nsys);
  if (docl), M.glob.dalpha = dV(end); end
    
end

   
%-------------------------------------------------------------------------------
function [Rcla, Ru_alpha, Rcla_U] = clalpha_residual(M)
% computes cl constraint (or just prescribed alpha) residual and Jacobian
% INPUT
%   M  : mfoil class with inviscid solution and post-processed cl_alpha, cl_ue
% OUTPUT
%   Rcla     : cl constraint residual = cl - cltgt (scalar)
%   Ru_alpha : lin of ue residual w.r.t. alpha (Nsys x 1)
%   Rcla_U   : lin of cl residual w.r.t state (1 x 4*Nsys)
% DETAILS
%   Used for cl-constrained mode, with alpha as the extra variable
%   Should be called with up-to-date cl and cl linearizations
  
  Nsys = M.glob.Nsys;   % number of dofs
  N = M.foil.N;         % number of points (dofs) on airfoil
  alpha = M.oper.alpha; % angle of attack (deg)
  
  if (M.oper.givencl)  % cl is prescribed, need to trim alpha    
    Rcla = M.post.cl - M.oper.cltgt; % cl constraint residual
    Rcla_U = [zeros(1,4*Nsys), M.post.cl_alpha];
    Rcla_U(4:4:4*N) = M.post.cl_ue; % only airfoil nodes affected 
    % Ru = ue - [uinv + ue_m*(ds.*ue)], uinv = uinvref*[cos(alpha);sin(alpha)]
    Ru_alpha = -get_ueinvref(M)*[-sind(alpha); cosd(alpha)]*pi/180;
  else     % alpha is prescribed, easy
    Rcla = 0; % no residual
    Ru_alpha = zeros(Nsys,1); % not really, but alpha is not changing
    Rcla_U = [zeros(1,4*Nsys), 1];
  end
  
end
  

%-------------------------------------------------------------------------------
function build_glob_sys(M)
% builds the primary variable global residual system for the coupled problem
% INPUT
%   M  : mfoil class with a valid solution in M.glob.U
% OUTPUT
%   M.glob.R   : global residual vector (3*Nsys x 1)
%   M.glob.R_U : residual Jacobian matrix (3*Nsys x 4*Nsys, sparse)
%   M.glob.R_x : residual linearization w.r.t. x (3*Nsys x Nsys, sparse)
% DETAILS
%   Loops over nodes/stations to assemble residual and Jacobian
%   Transition dicated by M.vsol.turb, which should be consistent with the state
%   Accounts for wake initialization and first-point similarity solutions
%   Also handles stagnation point on node via simple extrapolation
  
  Nsys = M.glob.Nsys;
  M.glob.R = zeros(3*Nsys,1);
  M.glob.R_U = sparse(3*Nsys,4*Nsys);
  M.glob.R_x = sparse(3*Nsys,Nsys);
  
  for is = 1:3  % loop over surfaces
    Is = M.vsol.Is{is}; % surface point indices
    xi = M.isol.xi(Is); % distance from LE stag point
    N = length(Is); % number of points on this surface
    U = M.glob.U(:,Is); % [th, ds, sa, ue] states at all points on this surface
    Aux = zeros(1,N); % auxiliary data at all points: [wgap]
    if (is<3),xg = M.foil.x(1,Is); end % global x locations of nodes on is
    xft = M.oper.xft(is)*M.geom.chord; % forced transition location
    
    % get parameter structure
    param = build_param(M, is);
    
    % set auxiliary data
    if (is == 3), Aux(1,:) = M.vsol.wgap; end
    
    % special case of tiny first xi -- will set to stagnation state later
    if (is < 3) && (xi(1) < 1e-8*xi(end)), i0 = 2;
    else, i0 = 1; % i0 indicates the "first" point station
    end
    
    % first point system 
    if (is < 3) 
      
      % calculate the stagnation state, a function of U1 and U2
      Ip = [i0,i0+1];
      [Ust, Ust_U, Ust_x, xst] = stagnation_state(U(:,Ip), xi(Ip)); % stag state
      param.turb = false; param.simi = true;  % similarity station flag  
      [R1, R1_Ut, ~] = residual_station(param, [xst,xst], [Ust, Ust], Aux(:,[i0,i0]));
      param.simi = false;
      R1_Ust = R1_Ut(:,1:4) + R1_Ut(:,5:8);
      R1_U = R1_Ust*Ust_U;
      R1_x = R1_Ust*Ust_x;
      J = [Is(i0), Is(i0+1)];
      
      if (i0 == 2) 
        % i0=1 point landed right on stagnation: set value to Ust
        vprint(param, 2, 'hit stagnation!\n');
        Ig = 3*Is(1) + (-2:0);
        Jg = 4*Is(1) + (-3:0);
        M.glob.R(Ig) = U(1:3,1) - Ust(1:3);
        M.glob.R_U(Ig,Jg) = M.glob.R_U(Ig,Jg) + eye(3,4);
        Jg = [4*J(1) + (-3:0), 4*J(2) + (-3:0)];
        M.glob.R_U(Ig,Jg) = M.glob.R_U(Ig,Jg) - Ust_U(1:3,:);
        M.glob.R_x(Ig,J) = -Ust_x(1:3,:);
      end
      
      % less accurate, first-order, original version
      %Ip = [i0,i0];
      %param.turb = false; param.simi = true;  % similarity station flag      
      %[R1, R1_U, R1_x] = residual_station(param, xi(Ip), U(:,Ip), Aux(:,Ip));
      %param.simi = false;
      %J = [Is(i0), Is(i0)]; % will add both contributions
      
    else
      % wake initialization
      [R1, R1_U, J] = wake_sys(M, param);
      R1_x = []; % no xi dependence of first wake residual
      param.turb = true; % force turbulent in wake if still laminar
      param.wake = true;
    end

    % store first point system in global residual, Jacobian
    Ig = 3*Is(i0) + (-2:0);
    M.glob.R(Ig) = R1;
    for j = 1:length(J)
      Jg = 4*J(j) + (-3:0);
      M.glob.R_U(Ig,Jg) = M.glob.R_U(Ig,Jg) + R1_U(:,4*j+(-3:0));
      if (~isempty(R1_x)), M.glob.R_x(Ig,J(j)) = M.glob.R_x(Ig,J(j)) + R1_x(:,j); end
    end
      
    % march over rest of points
    tran = false;
    for i = (i0+1):N
      Ip = [i-1,i]; % two points involved in the calculation
      
      % check for forced transition
      M.oper.forcet(is) = false;
      if (~tran) && (~param.turb) && (is<3) && (prod(xg(Ip)-xft)<0.)
        tran = true; 
        M.oper.xift(is) = xi(i-1) + (xi(i)-xi(i-1))*(xft-xg(i-1))/(xg(i)-xg(i-1));
        M.oper.forcet(is) = true;
      end;
      
      tran = xor(M.vsol.turb(Is(i-1)), M.vsol.turb(Is(i))); % transition flag
            
      % residual, Jacobian for point i
      if (tran)
        param.is = is;
        [Ri, Ri_U, Ri_x] = residual_transition(M, param, xi(Ip), U(:,Ip), Aux(:,Ip));
        store_transition(M, is, i);
      else
        [Ri, Ri_U, Ri_x] = residual_station(param, xi(Ip), U(:,Ip), Aux(:,Ip));
      end
        
      % store point i contribution in global residual, Jacobian
      Ig = 3*Is(i) + (-2:0); Jg = [4*Is(i-1)+(-3:0), 4*Is(i)+(-3:0)];
      M.glob.R(Ig) = M.glob.R(Ig) + Ri;
      M.glob.R_U(Ig,Jg) = M.glob.R_U(Ig,Jg) + Ri_U(:,1:8);
      M.glob.R_x(Ig,Is(Ip)) = M.glob.R_x(Ig,Is(Ip)) + Ri_x;
      
      % following transition, all stations will be turbulent
      if (tran), param.turb = true; end
      
    end
  end
  
  % special stagnation point treatment (not currently used)
  % residual_stagnation(M);
    
end


%-------------------------------------------------------------------------------
function [Ust, Ust_U, Ust_x, xst] = stagnation_state(U, x)
% extrapolates two states in U, first ones in BL, to stagnation
% INPUT
%   U  : [U1,U2] = states at first two nodes (4x2)
%   x  : [x1,x2] = x-locations of first two nodes (2x1)
% OUTPUT
%   Ust    : stagnation state (4x1)
%   Ust_U  : linearization of Ust w.r.t. U1 and U2 (4x8)
%   Ust_x  : linearization of Ust w.r.t. x1 and x2 (4x2)
%   xst    : stagnation point location ... close to 0
% DETAILS
%   fits a quadratic to the edge velocity: 0 at x=0, then through two states
%   linearly extrapolates other states in U to x=0, from U1 and U2

  % pull off states
  U1 = U(:,1); U2 = U(:,2); x1 = x(1); x2 = x(2);
  dx = x2-x1; dx_x = [-1, 1];
  rx = x2/x1; rx_x = [-rx,1]/x1;
  
  % linear extrapolation weights and stagnation state
  w1 =  x2/dx; w1_x = -w1/dx*dx_x + [ 0,1]/dx;
  w2 = -x1/dx; w2_x = -w2/dx*dx_x + [-1,0]/dx;
  Ust = U1*w1 + U2*w2;
  
  % quadratic extrapolation of the edge velocity for better slope, ue=K*x
  wk1 = rx/dx; wk1_x = rx_x/dx - wk1/dx*dx_x;
  wk2 = -1/(rx*dx); wk2_x = -wk2*(rx_x/rx + dx_x/dx); 
  K = wk1*U1(4) + wk2*U2(4);
  K_U = [0,0,0,wk1, 0,0,0,wk2];
  K_x = U1(4)*wk1_x + U2(4)*wk2_x;

  % less-accurate linear version
  %K = U1(4)/x1;
  %K_U = [0,0,0,1/x1, 0,0,0,0];
  %K_x = [-K/x1, 0];
  
  % stagnation coord cannot be zero, but must be small
  xst = 1e-6;
  Ust(4) = K*xst;  % linear dep of ue on x near stagnation
  Ust_U = [w1*eye(3,4), w2*eye(3,4); K_U*xst];
  Ust_x = [U1(1:3)*w1_x + U2(1:3)*w2_x; K_x*xst];
 
end


%-------------------------------------------------------------------------------
function [th, ds] = thwaites_init(K, nu)
% uses Thwaites correlation to initialize first node in stag point flow
% INPUT
%   K  : stagnation point constant
%   nu : kinematic viscosity
% OUTPUT
%   th : momentum thickness
%   ds : displacement thickness
% DETAILS
%   ue = K*x -> K = ue/x = stag point flow constant
%   th^2 = ue^(-6) * 0.45 * nu * int_0^x ue^5 dx = 0.45*nu/(6*K)
%   ds = Hstag*th = 2.2*th

  th = sqrt(0.45*nu/(6.*K)); % momentum thickness
  ds = 2.2*th; % displacement thickness

end


%-------------------------------------------------------------------------------
function [R, R_U, J] = wake_sys(M, param)
% constructs residual system corresponding to wake initialization
% INPUT
%   param  : parameters
% OUTPUT
%   R   : 3x1 residual vector for th, ds, sa
%   R_U : 3x12 residual linearization, as three 3x4 blocks
%   J   : indices of the blocks of U in R_U (lower, upper, wake)

  il = M.vsol.Is{1}(end); Ul = M.glob.U(:,il); % lower surface TE index, state
  iu = M.vsol.Is{2}(end); Uu = M.glob.U(:,iu); % upper surface TE index, state
  iw = M.vsol.Is{3}(  1); Uw = M.glob.U(:,iw); % first wake index, state
  [~, hTE, ~,~,~] = TE_info(M.foil.x); % trailing-edge gap is hTE

  % Obtain wake shear stress from upper/lower; transition if not turb
  param.turb = true; param.wake = false; % calculating turbulent quantities right before wake
  if (M.vsol.turb(il)), ctl = Ul(3); ctl_Ul = [0,0,1,0]; % already turb; use state
  else, [ctl, ctl_Ul] = get_cttr(Ul, param); end % transition shear stress, lower
  if (M.vsol.turb(iu)), ctu = Uu(3); ctu_Uu = [0,0,1,0]; % already turb; use state
  else, [ctu, ctu_Uu] = get_cttr(Uu, param); end % transition shear stress, upper
  thsum = Ul(1) + Uu(1); % sum of thetas
  ctw = (ctl*Ul(1) + ctu*Uu(1))/thsum; % theta-average
  ctw_Ul = (ctl_Ul*Ul(1) + (ctl - ctw)*[1,0,0,0])/thsum;
  ctw_Uu = (ctu_Uu*Uu(1) + (ctu - ctw)*[1,0,0,0])/thsum;

  % residual; note, delta star in wake includes the TE gap, hTE
  R = [Uw(1)-(Ul(1)+Uu(1)); Uw(2)-(Ul(2)+Uu(2)+hTE); Uw(3)-ctw];
  J = [il, iu, iw]; % R depends on states at these nodes
  R_Ul = [-eye(2,4); -ctw_Ul]; R_Uu = [-eye(2,4); -ctw_Uu]; R_Uw = eye(3,4);
  R_U = [R_Ul, R_Uu, R_Uw];

end
  

%-------------------------------------------------------------------------------
function Uw = wake_init(M, ue)
% initializes the first point of the wake, using data in M.glob.U
% INPUT
%   ue  : edge velocity at the wake point
% OUTPUT
%   Uw  : 4x1 state vector at the wake point
  
  iw = M.vsol.Is{3}(  1); Uw = M.glob.U(:,iw); % first wake index, state
  [R, ~, ~] = wake_sys(M, M.param); % construct the wake system
  Uw(1:3) = Uw(1:3) - R; Uw(4) = ue; % solve the wake system, use ue

end


%-------------------------------------------------------------------------------
function param = build_param(M, is)
% builds a parameter structure for side is
% INPUT
%   is  : side number, 1 = lower, 2 = upper, 3 = wake
% OUTPUT
%   param : M.param structure with side information

  param = M.param;  
  param.wake = (is == 3);
  param.turb = param.wake; % the wake is fully turbulent
  param.simi = false; % true for similarity station
  
end


%-------------------------------------------------------------------------------
function param = station_param(M, param, i)
% modifies parameter structure to be specific for station i
% INPUT
%   i  : station number (node index along the surface)
% OUTPUT
%   param : modified parameter structure
  param.turb = M.vsol.turb(i); % turbulent
  param.simi = ismember(i,M.isol.Istag); % similarity
end


%-------------------------------------------------------------------------------
function init_boundary_layer(M)
% initializes BL solution on foil and wake by marching with given edge vel, ue
% INPUT
%   The edge velocity field ue must be filled in on the airfoil and wake
% OUTPUT
%   The state in M.glob.U is filled in for each point

  Hmaxl = 3.8; % above this shape param value, laminar separation occurs
  Hmaxt = 2.5; % above this shape param value, turbulent separation occurs 
  
  ueinv = get_ueinv(M); % get inviscid velocity
  M.glob.Nsys = M.foil.N + M.wake.N; % number of global variables (nodes)
  
  % do we need to initialize?
  if (~M.oper.initbl) && (size(M.glob.U,2)==M.glob.Nsys)
    vprint(M.param,1, '\n <<< Starting with current boundary layer >>> \n');
    M.glob.U(4,:) = ueinv; % do set a new edge velocity
    return;
  end
    
  vprint(M.param,1, '\n <<< Initializing the boundary layer >>> \n');
  
  M.glob.U = zeros(4,M.glob.Nsys); % global solution matrix
  M.vsol.turb = zeros(M.glob.Nsys,1); % node flag: 0 = laminar, 1 = turbulent
  
  for is = 1:3  % loop over surfaces
    
    vprint(M.param, 3, '\nSide is = %d:\n', is);
    
    Is = M.vsol.Is{is}; % surface point indices
    xi = M.isol.xi(Is); % distance from LE stag point
    if (is<3),xg = M.foil.x(1,Is); end % global x locations of nodes on is
    ue = ueinv(Is); % edge velocities
    N = length(Is); % number of points
    U = zeros(4,N); % states at all points: [th, ds, sa, ue]
    Aux = zeros(1,N); % auxiliary data at all points: [wgap]
    xft = M.oper.xft(is)*M.geom.chord; % forced transition x location
    
    % ensure edge velocities are not tiny
    uemax = max(abs(ue)); ue = max(ue,1e-8*uemax);
     
    % get parameter structure
    param = build_param(M, is);
    
    % set auxiliary data
    if (is == 3), Aux(1,:) = M.vsol.wgap; end
    
    % initialize state at first point
    i0 = 1;
    if (is < 3) 

      % Solve for the stagnation state (Thwaites initialization + Newton)
      if (xi(1)<1e-8*xi(end)), K = ue(2)/xi(2); hitstag = true; 
      else, K = ue(1)/xi(1); hitstag = false;
      end
      [th, ds] = thwaites_init(K, param.mu0/param.rho0);
      xst = 1e-6; % small but nonzero
      Ust = [th; ds; 0; K*xst];
      nNewton = 20;
      for iNewton = 1:nNewton
        % call residual at stagnation
        param.turb = false; param.simi = true;  % similarity station flag 
        [R, R_U, ~] = residual_station(param, [xst,xst], [Ust,Ust], zeros(1,2));
        param.simi = false;
        if (norm(R) < 1e-10), break; end
        ID = 1:3; A = R_U(:, ID+4) + R_U(:,ID); b = -R; dU = [A\b; 0];
        % under-relaxation
        dm = max(abs([abs(dU(1)/Ust(1)), abs(dU(2)/Ust(2))]));
        omega = 1; if (dm > 0.2), omega = 0.2/dm; end
        dU = dU*omega;
        Ust = Ust + dU;
      end
      
      % store stagnation state in first one (rarely two) points
      if (hitstag), U(:,1) = Ust; U(4,1) = ue(1); i0=2; end
      U(:,i0) = Ust; U(4,i0) = ue(i0);
      
    else % wake
      U(:,1) = wake_init(M, ue(1)); % initialize wake state properly
      param.turb = true; % force turbulent in wake if still laminar
      M.vsol.turb(Is(1)) = true; % wake starts turbulent
    end
    
    % march over rest of points
    tran = false; % flag indicating that we are at transition
    i = i0+1;
    while (i<=N)
      Ip = [i-1,i]; % two points involved in the calculation
      
      % check for forced transition
      M.oper.forcet(is) = false;
      if (~tran) && (~param.turb) && (is<3) && (prod(xg(Ip)-xft)<0.)
        tran = true; 
        M.oper.xift(is) = xi(i-1) + (xi(i)-xi(i-1))*(xft-xg(i-1))/(xg(i)-xg(i-1));
        M.oper.forcet(is) = true;
        fprintf(1,'forced transition during marching: xft=%.5f, xift=%.5f\n', xft, M.oper.xift(is));
      end;
      U(:,i) = U(:,i-1); U(4,i) = ue(i); % guess = same state, new ue
      if (tran) % set shear stress at transition interval
        ct = get_cttr(U(:,i), param); U(3,i) = ct;
      end
      M.vsol.turb(Is(i)) = (tran || param.turb); % flag node i as turbulent
      direct = true; % default is direct mode
      nNewton = 30; iNswitch = 12;
      for iNewton = 1:nNewton
        
        % call residual at this station
        if (tran) % we are at transition
          vprint(param, 4, 'i=%d, residual_transition (iNewton = %d) \n', i, iNewton);
          try
            param.is = is;
            [R, R_U, ~] = residual_transition(M, param, xi(Ip), U(:,Ip), Aux(:,Ip));
          catch
            warning('Transition calculation failed in BL init. Continuing.');
            M.vsol.xt = 0.5*sum(xi(Ip));
            U(:,i) = U(:,i-1); U(4,i) = ue(i); U(3,i) = ct;
            R = 0; % so we move on
          end
        else
          vprint(param, 4, 'i=%d, residual_station (iNewton = %d)\n', i, iNewton);
          [R, R_U, ~] = residual_station(param, xi(Ip), U(:,Ip), Aux(:,Ip));
        end
        if (norm(R) < 1e-10), break; end
        
        if (direct) % direct mode => ue is prescribed => solve for th, ds, sa
          ID = 1:3; A = R_U(:, ID+4); b = -R; dU = [A\b; 0];
        else % inverse mode => Hk is prescribed 
          [Hk, Hk_U] = get_Hk(U(:,i), param);
          A = [R_U(:, 5:8); Hk_U]; b = [-R; Hktgt-Hk]; dU = A\b;
        end
          
        % under-relaxation
        dm = max(abs([abs(dU(1)/U(1,i-1)), abs(dU(2)/U(2,i-1))]));
        if (~direct), dm = max(dm, abs(dU(4)/U(4,i-1))); end
        if (param.turb), dm = max(dm, abs(dU(3)/U(3,i-1)));
        elseif (direct), dm = max(dm, abs(dU(3)/10)); 
        end
        omega = 1; if (dm > 0.3), omega = 0.3/dm; end
        dU = dU*omega;
        
        % trial update
        Ui = U(:,i) + dU;
        
        % clip extreme values
        if (param.turb), Ui(3) = max(min(Ui(3), .3), 1e-7); end
        %Hklim = 1.02; if (param.wake), Hklim = 1.00005; end
        %[Hk,Hk_U] = get_Hk(Ui, param);
        %dH = max(0,Hklim-Hk); Ui(2) = Ui(2) + dH*Ui(1);
        
        % check if about to separate
        Hmax = Hmaxl; if (param.turb), Hmax = Hmaxt; end
        [Hk,~] = get_Hk(Ui, param);

        if (direct) && ((Hk>Hmax) || (iNewton > iNswitch))
          % no update; need to switch to inverse mode: prescribe Hk
          direct = false;
          vprint(param, 2, '** switching to inverse: i=%d, iNewton=%d\n', i, iNewton);
          [Hk,~] = get_Hk(U(:,i-1), param); Hkr = (xi(i)-xi(i-1))/U(1,i-1);
          if (param.wake)
            H2 = Hk; 
            for k=1:6, H2 = H2 - (H2+.03*Hkr*(H2-1)^3-Hk)/(1+.09*Hkr*(H2-1)^2); end
            Hktgt = max(H2, 1.01);
          elseif (param.turb), Hktgt = Hk - .15*Hkr; % turb: decrease in Hk
          else, Hktgt = Hk + .03*Hkr; % lam: increase in Hk 
          end
          if (~param.wake), Hktgt = max(Hktgt, Hmax); end
          if (iNewton > iNswitch), U(:,i) = U(:,i-1); U(4,i) = ue(i); end % reinit
        else, U(:,i) = Ui;  % take the update
        end
      end
      if (iNewton >= nNewton)
        vprint(param, 1, '** BL init not converged: is=%d, i=%d **\n', is, i); 
        % extrapolate values
        U(:,i) = U(:,i-1); U(4,i) = ue(i);
        if (is<3)
          U(1,i) = U(1,i-1)*(xi(i)/xi(i-1))^.5;
          U(2,i) = U(2,i-1)*(xi(i)/xi(i-1))^.5;
        else
          rlen = (xi(i)-xi(i-1))/(10*U(2,i-1));
          U(2,i) = (U(2,i-1) + U(1,i-1)*rlen)/(1+rlen);
        end
      end
      
      % check for transition
      if (~param.turb) && (~tran) && (U(3,i)>param.ncrit)
        vprint(param,2, 'Identified transition at (is=%d, i=%d): n=%.5f, ncrit=%.5f\n', ...
               is, i, U(3,i), param.ncrit);
        tran = true; 
        continue; % redo station with transition
      end
      
      if (tran)
        store_transition(M, is, i);  % store transition location
        param.turb = true; tran = false; % turbulent after transition
        fprintf(1, 'storing transitinon\n');
      end 
      
      i = i+1; % next point
    end
    
    % store states
    M.glob.U(:,Is) = U;
    
  end
  
end


%-------------------------------------------------------------------------------
function store_transition(M, is, i)
% stores xi and x transition locations using current M.vsol.xt 
% INPUT
%   is,i : side,station number
% OUTPUT
%   M.vsol.Xt stores the transition location s and x values
  
  forcet = M.oper.forcet(is); % flag to force transition
  xft = M.oper.xft(is)*M.geom.chord; % forced transition x location
  i0 = M.vsol.Is{is}(i-1); i1 = M.vsol.Is{is}(i); % pre/post transition nodes
  assert((i0<=M.foil.N) && (i1<=M.foil.N), 'Can only store transition on airfoil');
  xi0 = M.isol.xi(i0); xi1 = M.isol.xi(i1); % xi (s) locations at nodes
  x0 = M.foil.x(1,i0); x1 = M.foil.x(1,i1); % x locations at nodes

  xift = xi0 + (xi1-xi0)*(xft-x0)/(x1-x0);
  if ((~forcet) || ((M.vsol.xt > 0) && (M.vsol.xt < xift)))
    xt = M.vsol.xt; spre = 'free'; % free transition
  else
    xt = xift; M.oper.xift(is) = xt; spre = 'forced';
  end
  if ((xt<xi0) || (xt>xi1))
    vprint(M.param,1, 'Warning: transition (%.3f) off interval (%.3f,%.3f)!\n', xt, xi0, xi1);
  end
  M.vsol.Xt(is,1) = xt; % xi location
  M.vsol.Xt(is,2) = x0 + (xt-xi0)/(xi1-xi0)*(x1-x0); % x location
  slu = {'lower', 'upper'};
  vprint(M.param,1, '  %s transition on %s side at x=%.5f\n', spre, slu{is}, M.vsol.Xt(is,2));

end
   

%-------------------------------------------------------------------------------
function update_transition(M)
% updates transition location using current state
% INPUT
%   a valid state in M.glob.U
% OUTPUT
%   M.vsol.turb : updated with latest lam/turb flags for each node
%   M.glob.U    : updated with amp factor or shear stress as needed at each node
  
  for is = 1:2  % loop over lower/upper surfaces
    
    Is = M.vsol.Is{is}; % surface point indices
    N = length(Is); % number of points
    xft = M.oper.xft(is)*M.geom.chord; % forced transition x location

    % get parameter structure
    param = build_param(M, is);
    
    % current last laminar station
    I = find(M.vsol.turb(Is) == 0); ilam0 = I(end);
    
    % current amp/ctau solution (so we do not change it unnecessarily)
    sa = M.glob.U(3,Is);
    
    % march amplification equation to get new last laminar station
    ilam = march_amplification(M, is);
    
    % set xi value
    if (M.oper.forcet(is)), 
      i0 = Is(ilam); i1 = Is(ilam+1);
      M.oper.xift(is) = M.isol.xi(i0) + (M.isol.xi(i1)-M.isol.xi(i0))*(xft-M.foil.x(1,i0))/(M.foil.x(1,i1)-M.foil.x(1,i0));
    end
    
    if (ilam == ilam0), M.glob.U(3,Is) = sa; continue; end % no change
    
    vprint(param, 2, '  Update transition: last lam [%d]->[%d]\n', ilam0, ilam);
    
    if (ilam < ilam0)
      % transition is now earlier: fill in turb between [ilam+1, ilam0]
      param.turb = true;
      [sa0, ~] = get_cttr(M.glob.U(:,Is(ilam+1)), param);
      sa1 = sa0; if (ilam0<N), sa1 = M.glob.U(3,Is(ilam0+1)); end
      xi = M.isol.xi(Is); dx = xi(min(ilam0+1,N))-xi(ilam+1);
      for i = (ilam+1):ilam0
        if (dx==0) || (i==ilam+1), f = 0; else, f = (xi(i)-xi(ilam+1))/dx; end
        if ((ilam+1) == ilam0), f = 1; end
        M.glob.U(3,Is(i)) = sa0 + f*(sa1-sa0);
        assert(M.glob.U(3,Is(i)) > 0, 'negative ctau in update_transition');
        M.vsol.turb(Is(i)) = 1;
      end      

    elseif (ilam > ilam0)
      % transition is now later: lam already filled in; leave turb alone
      for i = ilam0:ilam, M.vsol.turb(Is(i)) = 0; end
    end
  end 
  
end


%-------------------------------------------------------------------------------
function [ilam] = march_amplification(M, is)
% marches amplification equation on surface is
% INPUT
%   is : surface number index
% OUTPUT
%   ilam : index of last laminar station before transition
%   M.glob.U : updated with amp factor at each (new) laminar station
  
  Is = M.vsol.Is{is}; % surface point indices
  N = length(Is); % number of points
  param = build_param(M, is); % get parameter structure
  U = M.glob.U(:,Is); % states
  turb = M.vsol.turb(Is); % turbulent station flag
  xft = M.oper.xft(is)*M.geom.chord; % forced transition x location
  
  % loop over stations, calculate amplification
  U(3,1) = 0.; % no amplification at first station
  param.turb = false; param.wake = false;
  i = 2;
  while (i <= N)
    U1 = U(:,i-1); U2 = U(:,i); % states
    if (turb(i)), U2(3) = U1(3)*1.01; end % initialize amp if turb
    dx = M.isol.xi(Is(i)) - M.isol.xi(Is(i-1)); % interval length
        
    % Newton iterations, only needed if adding extra amplification in damp
    nNewton = 20;
    for iNewton = 1:nNewton
      % amplification rate, averaged
      [damp1, damp1_U1] = get_damp(U1, param);
      [damp2, damp2_U2] = get_damp(U2, param);
      [damp, damp_U] = upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2);
    
      Ramp = U2(3) - U1(3) - damp*dx;
      
      if (iNewton > 12)
        vprint(param,3,'i=%d, iNewton=%d, sa = [%.5e, %.5e], damp = %.5e, Ramp = %.5e\n', ...
                i, iNewton, U1(3), U2(3), damp, Ramp);
      end
      
      if (abs(Ramp)<1e-12), break; end % converged
      Ramp_U = [0,0,-1,0, 0,0,1,0] - damp_U*dx;
      dU = -Ramp/Ramp_U(7);
      omega = 1; dmax = 0.5*(1.01-iNewton/nNewton);
      if (abs(dU) > dmax), omega = dmax/abs(dU); end
      U2(3) = U2(3) + omega*dU;
    end
    if (iNewton >= nNewton), vprint(param,1, 'march amp Newton unconverged!\n'); end
    
    % check for transition
    M.oper.forcet(is) = false;
    if (prod(M.foil.x(1,Is([i-1,i]))-xft)<0.) 
      M.oper.forcet(is) = true;
      vprint(param, 2,'  forced transition (is,i=%d,%d)\n', is, i);
      break
    elseif (U2(3)>param.ncrit)
      vprint(param, 2,'  march_amplification (is,i=%d,%d): %.5e is above critical.\n', ...
             is, i, U2(3));
      break
    else
      M.glob.U(3,Is(i)) = U2(3); % store amplification in M.glob.U      
      U(3,i) = U2(3); % also store in local copy!
      if (imag(U(3,i))), error('imaginary amp during march'); end
    end
    
    i = i+1; % next station
  end
  
  ilam = i-1; % set last laminar station
  
end


%-------------------------------------------------------------------------------
function [R, R_U, R_x] = residual_transition(M, param, x, U, Aux)
% calculates the combined lam + turb residual for a transition station
% INPUT
%   param : parameter structure
%   x     : 2x1 vector, [x1, x2], containing xi values at the points
%   U     : 4x2 matrix, [U1, U2], containing the states at the points
%   Aux   : ()x2 matrix, [Aux1, Aux2] of auxiliary data at the points
% OUTPUT
%   R     : 3x1 transition residual vector
%   R_U   : 3x8 residual Jacobian, [R_U1, R_U2]
%   R_x   : 3x2 residual linearization w.r.t. x, [R_x1, R_x2]
% DETAILS
%   The state U1 should be laminar; U2 should be turbulent
%   Calculates and linearizes the transition location in the process
%   Assumes linear variation of th and ds from U1 to U2
  
  % states
  U1 = U(:,1); U2 = U(:,2); sa = U(3,:);
  I1 = 1:4; I2 = 5:8; Z = zeros(1,4);

  % forced transition?
  is = param.is;
  forcet = M.oper.forcet(is);  % flag to force transition
  xft = M.oper.xft(is)*M.geom.chord; % forced transition x location
  xift = M.oper.xift(is); % forced transition xi location
  
  % interval
  x1 = x(1); x2 = x(2); dx = x2-x1;
  
  if (~forcet)
    % determine transition location (xt) using amplification equation
    xt = x1 + 0.5*dx; % guess
    ncrit = param.ncrit; % critical amp factor
    nNewton = 20;
    vprint(param, 3, '  Transition interval = [%.5e, %.5e]\n', x1, x2);
    %  U1, U2
    for iNewton = 1:nNewton
      w2 = (xt-x1)/dx; w1 = 1-w2; % weights
      Ut = w1*U1 + w2*U2; Ut_xt = (U2-U1)/dx; % state at xt
      Ut(3) = ncrit; Ut_xt(3) = 0.; % amplification at transition
      [damp1, damp1_U1] = get_damp(U1, param);
      [dampt, dampt_Ut] = get_damp(Ut, param); dampt_Ut(3) = 0.;
      Rxt = ncrit - sa(1) - 0.5*(xt-x1)*(damp1 + dampt);
      Rxt_xt = -0.5*(damp1+dampt) - 0.5*(xt-x1)*(dampt_Ut*Ut_xt);
      dxt = -Rxt/Rxt_xt;
      vprint(param, 4, '   Transition: iNewton,Rxt,xt = %d,%.5e,%.5e\n', iNewton,Rxt,xt);
      dmax = 0.2*dx*(1.1-iNewton/nNewton);
      if (abs(dxt)>dmax), dxt = dxt*dmax/abs(dxt); end
      if (abs(Rxt) < 1e-10), break; end
      if (iNewton<nNewton), xt = xt + dxt; end
    end

    if (iNewton >= nNewton), warning('Transition location calculation failed.'); end

    % prepare for xt linearizations
    Rxt_U = -0.5*(xt-x1)*[damp1_U1 + dampt_Ut*w1, dampt_Ut*w2]; Rxt_U(3) = Rxt_U(3)-1;
    Ut_x1 = (U2-U1)*(w2-1)/dx; Ut_x2 = (U2-U1)*(-w2)/dx; % at fixed xt
    Ut_x1(3) = 0; Ut_x2(3) = 0; % amp at xt is always ncrit
    Rxt_x1 = 0.5*(damp1+dampt) - 0.5*(xt-x1)*(dampt_Ut*Ut_x1);
    Rxt_x2 =                   - 0.5*(xt-x1)*(dampt_Ut*Ut_x2);
    
    % sensitivity of xt w.r.t. U,x from Rxt(xt,U,x) = 0 constraint
    xt_U = -Rxt_U/Rxt_xt;  xt_U1 = xt_U(:,I1); xt_U2 = xt_U(:,I2);
    xt_x1 = -Rxt_x1/Rxt_xt; xt_x2 = -Rxt_x2/Rxt_xt;
  
  else
    xt = xift;  % transition xi location
    w2 = (xt-x1)/dx; w1 = 1-w2; % weights
    Ut = w1*U1 + w2*U2; Ut_xt = (U2-U1)/dx; % state at xt
    Rxt = 0; Rxt_xt = 1; Rxt_x1 = -w1; Rxt_x2 = -w2;
    Rxt_U = zeros(1,8);
    Ut_x1 = (U2-U1)*(w2-1)/dx; Ut_x2 = (U2-U1)*(-w2)/dx; % at fixed xt
    xt_U = -Rxt_U/Rxt_xt;  xt_U1 = xt_U(:,I1); xt_U2 = xt_U(:,I2);
    xt_x1 = -Rxt_x1/Rxt_xt; xt_x2 = -Rxt_x2/Rxt_xt;
    
  end
  M.vsol.xt = xt; % save transition location
  
  % include derivatives w.r.t. xt in Ut_x1 and Ut_x2
  Ut_x1 = Ut_x1 + Ut_xt*xt_x1;
  Ut_x2 = Ut_x2 + Ut_xt*xt_x2;
  
  % sensitivity of Ut w.r.t. U1 and U2
  Ut_U1 = w1*eye(4) + (U2-U1)*xt_U1/dx; % w1*I + U1*w1_xt*xt_U1 + U2*w2_xt*xt_U1;
  Ut_U2 = w2*eye(4) + (U2-U1)*xt_U2/dx; % w2*I + U1*w1_xt*xt_U2 + U2*w2_xt*xt_U2;
  
  % laminar and turbulent states at transition
  Utl = Ut; Utl_U1 = Ut_U1; Utl_U2 = Ut_U2; Utl_x1 = Ut_x1; Utl_x2 = Ut_x2;
  if (~forcet), Utl(3) = ncrit; Utl_U1(3,:) = Z; Utl_U2(3,:) = Z; Utl_x1(3) = 0; Utl_x2(3) = 0; end;
  Utt = Ut; Utt_U1 = Ut_U1; Utt_U2 = Ut_U2; Utt_x1 = Ut_x1; Utt_x2 = Ut_x2;
  
  % parameter structure
  param = build_param(M, 0);
  
  % set turbulent shear coefficient, sa, in Utt
  param.turb = true;
  [cttr, cttr_Ut] = get_cttr(Ut, param);
  Utt(3) = cttr; Utt_U1(3,:) = cttr_Ut*Ut_U1; Utt_U2(3,:) = cttr_Ut*Ut_U2;
  Utt_x1(3) = cttr_Ut*Ut_x1; Utt_x2(3) = cttr_Ut*Ut_x2;
  
  % laminar/turbulent residuals and linearizations
  param.turb = false;
  [Rl, Rl_U, Rl_x] = residual_station(param, [x1,xt], [U1,Utl], Aux);
  Rl_U1 = Rl_U(:,I1); Rl_Utl = Rl_U(:,I2);
  if (forcet), Rl(3) = 0; Rl_U1(3,:) = Z; Rl_Utl(3,:) = Z; end; % !!
  param.turb = true;
  [Rt, Rt_U, Rt_x] = residual_station(param, [xt,x2], [Utt,U2], Aux);
  Rt_Utt = Rt_U(:,I1); Rt_U2 = Rt_U(:,I2);

  % combined residual and linearization
  R = Rl + Rt;
  if (any(imag(R))), error('imaginary transition residual');end
  R_U1 = Rl_U1 + Rl_Utl*Utl_U1 + Rl_x(:,2)*xt_U1 + Rt_Utt*Utt_U1 + Rt_x(:,1)*xt_U1;
  R_U2 = Rl_Utl*Utl_U2 + Rl_x(:,2)*xt_U2 + Rt_Utt*Utt_U2 + Rt_U2 + Rt_x(:,1)*xt_U2;
  R_U = [R_U1, R_U2];
  R_x = [Rl_x(:,1) + Rl_x(:,2)*xt_x1 + Rt_x(:,1)*xt_x1 + Rl_Utl*Utl_x1 + Rt_Utt*Utt_x1, ...
         Rt_x(:,2) + Rl_x(:,2)*xt_x2 + Rt_x(:,1)*xt_x2 + Rl_Utl*Utl_x2 + Rt_Utt*Utt_x2];
  
end


%-------------------------------------------------------------------------------
function [R, R_U, R_x] = residual_station(param, x, U, Aux)
% calculates the viscous residual at one non-transition station
% INPUT
%   param : parameter structure
%   x     : 2x1 vector, [x1, x2], containing xi values at the points
%   U     : 4x2 matrix, [U1, U2], containing the states at the points
%   Aux   : ()x2 matrix, [Aux1, Aux2] of auxiliary data at the points
% OUTPUT
%   R     : 3x1 residual vector (mom, shape-param, amp/lag)
%   R_U   : 3x8 residual Jacobian, [R_U1, R_U2]
%   R_x   : 3x2 residual linearization w.r.t. x, [R_x1, R_x2]
% DETAILS
%   The input states are U = [U1, U2], each with th,ds,sa,ue

  % modify ds to take out wake gap (in Aux) for all calculations below
  U(2,:) = U(2,:) - Aux(1,:);
  
  % states
  U1 = U(:,1); U2 = U(:,2); Um = 0.5*(U1+U2);
  th = U(1,:); ds = U(2,:); sa = U(3,:); 
  
  % speed needs compressibility correction
  [uk1, uk1_u] = get_uk(U1(4),param);
  [uk2, uk2_u] = get_uk(U2(4),param);
    
  % log changes
  thlog = log(th(2)/th(1));
  thlog_U = [-1/th(1),0,0,0, 1/th(2),0,0,0];
  uelog = log(uk2/uk1);
  uelog_U = [0,0,0,-uk1_u/uk1, 0,0,0,uk2_u/uk2];
  xlog = log(x(2)/x(1)); xlog_x = [-1/x(1), 1/x(2)];
  dx = x(2)-x(1); dx_x = [-1, 1];
  
  % upwinding factor
  [upw, upw_U] = get_upw(U1, U2, param);
  
  % shape parameter
  [H1, H1_U1] = get_H(U(:,1));
  [H2, H2_U2] = get_H(U(:,2));
  H = 0.5*(H1+H2);
  H_U = 0.5*[H1_U1, H2_U2];
  
  % Hstar = KE shape parameter, averaged
  [Hs1, Hs1_U1] = get_Hs(U1, param);
  [Hs2, Hs2_U2] = get_Hs(U2, param);
  [Hs, Hs_U] = upwind(0.5, 0, Hs1, Hs1_U1, Hs2, Hs2_U2);  
  
  % log change in Hstar
  Hslog = log(Hs2/Hs1);
  Hslog_U = [-1/Hs1*Hs1_U1, 1/Hs2*Hs2_U2];
  
  % similarity station is special: U1 = U2, x1 = x2
  if (param.simi)
    thlog = 0; thlog_U = thlog_U*0;
    Hslog = 0; Hslog_U = Hslog_U*0;
    uelog = 1; uelog_U = uelog_U*0;
    xlog = 1; xlog_x = [0, 0];
    dx = 0.5*(x(1)+x(2)); dx_x = [0.5,0.5];
  end
    
  % Hw = wake shape parameter
  [Hw1, Hw1_U1] = get_Hw(U(:,1), Aux(1,1));
  [Hw2, Hw2_U2] = get_Hw(U(:,2), Aux(1,2));
  Hw = 0.5*(Hw1 + Hw2);
  Hw_U = 0.5*[Hw1_U1, Hw2_U2];
  
  % set up shear lag or amplification factor equation
  if (param.turb)

    % log change of root shear stress coeff
    salog = log(sa(2)/sa(1));
    salog_U = [0,0,-1/sa(1),0, 0,0,1/sa(2),0];
    
    % BL thickness measure, averaged
    [de1, de1_U1] = get_de(U1, param);
    [de2, de2_U2] = get_de(U2, param);
    [de, de_U] = upwind(0.5, 0, de1, de1_U1, de2, de2_U2);
    
    % normalized slip velocity, averaged
    [Us1, Us1_U1] = get_Us(U1, param);
    [Us2, Us2_U2] = get_Us(U2, param);
    [Us, Us_U] = upwind(0.5, 0, Us1, Us1_U1, Us2, Us2_U2);
    
    % Hk, upwinded
    [Hk1, Hk1_U1] = get_Hk(U1, param);
    [Hk2, Hk2_U2] = get_Hk(U2, param);
    [Hk, Hk_U] = upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2);
    
    % Re_theta, averaged
    [Ret1, Ret1_U1] = get_Ret(U1, param);
    [Ret2, Ret2_U2] = get_Ret(U2, param);
    [Ret, Ret_U] = upwind(0.5, 0, Ret1, Ret1_U1, Ret2, Ret2_U2);
    
    % skin friction, upwinded
    [cf1, cf1_U1] = get_cf(U1, param);
    [cf2, cf2_U2] = get_cf(U2, param);
    [cf, cf_U] = upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2);
    
    % displacement thickness, averaged
    dsa = 0.5*(ds(1) + ds(2));
    dsa_U = 0.5*[0,1,0,0, 0,1,0,0];
    
    % uq = equilibrium 1/ue * due/dx
    [uq, uq_U] = get_uq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param);
    
    % cteq = root equilibrium wake layer shear coeficient: (ctau eq)^.5
    [cteq1, cteq1_U1] = get_cteq(U1, param);
    [cteq2, cteq2_U2] = get_cteq(U2, param);
    [cteq, cteq_U] = upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2);

    % root of shear coefficient (a state), upwinded
    [saa, saa_U] = upwind(upw, upw_U, sa(1), [0,0,1,0], sa(2), [0,0,1,0]);
    
    % lag coefficient
    Klag = param.SlagK;
    beta = param.GB;
    Clag = Klag/beta*1/(1+Us);
    Clag_U = -Clag/(1+Us)*Us_U;
    
    % extra dissipation in wake
    ald = 1.0;
    if (param.wake), ald = param.Dlr; end
    
    % shear lag equation
    Rlag = Clag*(cteq-ald*saa)*dx - 2*de*salog + 2*de*(uq*dx-uelog)*param.Cuq;
    Rlag_U = Clag_U*(cteq-ald*saa)*dx + Clag*(cteq_U-ald*saa_U)*dx ...
             - 2*de_U*salog - 2*de*salog_U ...
             + 2*de_U*(uq*dx-uelog)*param.Cuq + 2*de*(uq_U*dx-uelog_U)*param.Cuq;
    Rlag_x = Clag*(cteq-ald*saa)*dx_x + 2*de*uq*dx_x;
        
  else
    % laminar, amplification factor equation
    
    if (param.simi)
      % similarity station
      Rlag = sa(1) + sa(2); % no amplification
      Rlag_U = [0,0,1,0, 0,0,1,0];
      Rlag_x = [0,0];
    else
      % amplification factor equation in Rlag

      % amplification rate, averaged
      [damp1, damp1_U1] = get_damp(U1, param);
      [damp2, damp2_U2] = get_damp(U2, param);
      [damp, damp_U] = upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2);
      
      Rlag = sa(2) - sa(1) - damp*dx;
      Rlag_U = [0,0,-1,0, 0,0,1,0] - damp_U*dx;
      Rlag_x = -damp*dx_x;
    end
    
  end 
  
  % squared mach number, symmetrical average
  [Ms1, Ms1_U1] = get_Mach2(U1, param);
  [Ms2, Ms2_U2] = get_Mach2(U2, param);
  [Ms, Ms_U] = upwind(0.5, 0, Ms1, Ms1_U1, Ms2, Ms2_U2);
  
  % skin friction * x/theta, symmetrical average
  [cfxt1, cfxt1_U1, cfxt1_x1] = get_cfxt(U1, x(1), param);
  [cfxt2, cfxt2_U2, cfxt2_x2] = get_cfxt(U2, x(2), param);
  [cfxtm, cfxtm_Um, cfxtm_xm] = get_cfxt(Um, 0.5*(x(1)+x(2)), param);
  cfxt = 0.25*cfxt1 + 0.5*cfxtm + 0.25*cfxt2;
  cfxt_U = 0.25*[cfxt1_U1+cfxtm_Um, cfxtm_Um+cfxt2_U2];
  cfxt_x = 0.25*[cfxt1_x1+cfxtm_xm, cfxtm_xm+cfxt2_x2];
  
  % momentum equation
  Rmom = thlog + (2+H+Hw-Ms)*uelog - 0.5*xlog*cfxt;
  Rmom_U = thlog_U + (H_U+Hw_U-Ms_U)*uelog + (2+H+Hw-Ms)*uelog_U - 0.5*xlog*cfxt_U;
  Rmom_x = -0.5*xlog_x*cfxt - 0.5*xlog*cfxt_x;  

  % dissipation function times x/theta: cDi = (2*cD/H*)*x/theta, upwinded
  [cDixt1, cDixt1_U1, cDixt1_x1] = get_cDixt(U1, x(1), param);
  [cDixt2, cDixt2_U2, cDixt2_x2] = get_cDixt(U2, x(2), param);
  [cDixt, cDixt_U] = upwind(upw, upw_U, cDixt1, cDixt1_U1, cDixt2, cDixt2_U2);
  cDixt_x = [(1-upw)*cDixt1_x1, upw*cDixt2_x2];
  
  % cf*x/theta, upwinded
  [cfxtu, cfxtu_U] = upwind(upw, upw_U, cfxt1, cfxt1_U1, cfxt2, cfxt2_U2);
  cfxtu_x = [(1-upw)*cfxt1_x1, upw*cfxt2_x2];
  
  % Hss = density shape parameter, averaged
  [Hss1, Hss1_U1] = get_Hss(U1, param);
  [Hss2, Hss2_U2] = get_Hss(U2, param);
  [Hss, Hss_U] = upwind(0.5, 0, Hss1, Hss1_U1, Hss2, Hss2_U2);
  
  Rshape = Hslog + (2*Hss/Hs + 1-H-Hw)*uelog + xlog*(0.5*cfxtu - cDixt);
  Rshape_U = Hslog_U + (2*Hss_U/Hs - 2*Hss/Hs^2*Hs_U -H_U - Hw_U)*uelog + ...
      (2*Hss/Hs + 1-H-Hw)*uelog_U + xlog*(0.5*cfxtu_U - cDixt_U);
  Rshape_x = xlog_x*(0.5*cfxtu - cDixt) + xlog*(0.5*cfxtu_x - cDixt_x);
  
  % put everything together
  R = [Rmom; Rshape; Rlag];
  R_U = [Rmom_U; Rshape_U; Rlag_U];
  R_x = [Rmom_x; Rshape_x; Rlag_x];

end


%-------------------------------------------------------------------------------
function residual_stagnation(M)
% replaces the residual and Jacobian in the global system with the 
% stagnation residual and Jacobian, at the two stagnation points
 
  param = build_param(M, 1); % parameters
  param.turb = false; param.simi = false; param.wake = false;

  
  Ist = M.isol.Istag; % stagnation point indices
  I = [Ist(1)-1, Ist(1), Ist(2), Ist(2)+1]; % surrounding points too
  x = M.isol.xi(I); % x-coords, measured away from stag
  U = M.glob.U(:,I); % 4 states

  % weights for calculating the stagnation state
  dx = x(3)+x(2); dx_x = [1,1]; % derivatives refer to points 2,3
  w2 =  x(3)/dx; w2_x = -w2/dx*dx_x + [ 0,1]/dx;
  w3 =  x(2)/dx; w3_x = -w3/dx*dx_x + [ 1,0]/dx;
  
  % stagnation state
  Ust = w2*U(:,2) + w3*U(:,3);
  K = 0.5*(U(4,2)/x(2) + U(4,3)/x(3)); % ue = K*x at stag
  K_U = [0,0,0,0.5/x(2), 0,0,0,0.5/x(3)];
  K_x = 0.5*[-U(4,2)/x(2)^2, -U(4,3)/x(3)^2];
  xst = 1e-6*M.geom.chord; % xst needs to be small but nonzero
  Ust(4) = K*xst;
  Ust_U = [w2*eye(3,4), w3*eye(3,4); K_U*xst];
  Ust_x = [U(1:3,2)*w2_x + U(1:3,3)*w3_x; K_x*xst];
  
  % ue and x log quantities at stagnation (both have to be the same constant)
  uelog = 1;
  xlog = 1;
  
  % shape parameter
  [H, H_Ust] = get_H(Ust);
  
  % squared Mach number
  [Ms, Ms_Ust] = get_Mach2(Ust, param);
  
  % skin friction * x/theta
  [cfxt, cfxt_Ust, ~] = get_cfxt(Ust, xst, param);
  
  % momentum equation at stagnation
  Rmom = (2+H-Ms)*uelog - 0.5*xlog*cfxt;
  Rmom_Ust = (H_Ust-Ms_Ust)*uelog - 0.5*xlog*cfxt_Ust;
  Rmom_U = Rmom_Ust*Ust_U;
  Rmom_x = Rmom_Ust*Ust_x;
  
  % dissipation function times x/theta: cDi = (2*cD/H*)*x/theta
  [cDixt, cDixt_Ust, ~] = get_cDixt(Ust, xst, param);
    
  % Hstar = KE shape parameter, averaged
  [Hs, Hs_Ust] = get_Hs(Ust, param);
  
  % Hss = density shape parameter
  [Hss, Hss_Ust] = get_Hss(Ust, param);
  
  % shape parameter equation at stagnation
  Rshape = (2*Hss/Hs+1-H)*uelog + xlog*(0.5*cfxt - cDixt);
  Rshape_Ust = (2*Hss_Ust/Hs - 2*Hss/Hs^2*Hs_Ust - H_Ust)*uelog + ...
      xlog*(0.5*cfxt_Ust - cDixt_Ust);
  Rshape_U = Rshape_Ust*Ust_U;
  Rshape_x = Rshape_Ust*Ust_x;
  
  % amplification equation at stagnation
  Ramp = Ust(3); % no amplification
  Ramp_Ust = [0,0,1,0];
  Ramp_U = Ramp_Ust*Ust_U;
  Ramp_x = Ramp_Ust*Ust_x;
  
  % put stagnation residual into the global system -- Ist(1) eqn
  Ig = 3*Ist(1) + (-2:0); Jg = [4*I(2)+(-3:0), 4*I(3)+(-3:0)];
  M.glob.R(Ig) = [Rmom; Rshape; Ramp];
  M.glob.R_U(Ig,:) = 0; M.glob.R_x(Ig,:) = 0; % clear out rows
  M.glob.R_U(Ig,Jg) = [Rmom_U; Rshape_U; Ramp_U];
  M.glob.R_x(Ig,Ist) = [Rmom_x; Rshape_x; Ramp_x];
  
  % second equation: second-order BL eqns between points 2 and 3

  % CKx = C/K*dx; ue = K*x + 0.5*C*x^2
  ue = U(4,:);
  u12 = ue(1)/ue(2);
  u43 = ue(4)/ue(3);
  x21 = x(2)/x(1); x21_x = [-x21,1,0,0]/x(1);
  x34 = x(3)/x(4); x34_x = [0,0,1,-x34]/x(4);
  dx = x(3)+x(2); dx_x = [0,1,1,0]; % stag interval length; note x1,x2 are neg
  dx1 = dx/x(1); dx1_x = dx_x/x(1)-dx/x(1)^2*[1,0,0,0];
  dx4 = dx/x(4); dx4_x = dx_x/x(4)-dx/x(4)^2*[0,0,0,1];
  CKx = (1-u12*x21)*dx1 - (1-u43*x34)*dx4;
  CKx_ue = [ [-1, u12]/ue(2)*x21*dx1, [-u43,1]/ue(3)*x34*dx4 ];
  CKx_x = (-u12*x21_x)*dx1 + (1-u12*x21)*dx1_x ...
          - (-u43*x34_x)*dx4 - (1-u43*x34)*dx4_x;
    
  % momentum equation
  m0 = (2+H-Ms); % at stagnation
  m0_Ust = (H_Ust-Ms_Ust);
  th0 = Ust(1); th0_Ust = [1,0,0,0];
  dth = U(1,3)-U(1,2); dth_U = [-1,0,0,0, 1,0,0,0];
  [H2, H2_U2] = get_H(U(:,2));
  [H3, H3_U3] = get_H(U(:,3));
  [Ms2, Ms2_U2] = get_Mach2(U(:,2), param);
  [Ms3, Ms3_U3] = get_Mach2(U(:,3), param);
  dm = (2+H3-Ms3)-(2+H2-Ms2); dm_U = [-H2_U2+Ms2_U2, H3_U3-Ms3_U3];
  [F2, F2_U2] = get_cfutstag(U(:,2), param);
  [F3, F3_U3] = get_cfutstag(U(:,3), param);
  dF = F3-F2; dF_U = [-F2_U2, F3_U3];
  Rmom = (1+2*m0)*dth/th0 + dm + m0*CKx - 0.5*dF/(K*th0^2);
  Rmom_Ust = 2*m0_Ust*dth/th0 - (1+2*m0)*dth/th0^2*th0_Ust ...
      + m0_Ust*CKx + dF/(K*th0^3)*th0_Ust;
  Rmom_U = (1+2*m0)*dth_U/th0 + dm_U - 0.5*dF_U/(K*th0^2) ...
           + 0.5*dF/(K^2*th0^2)*K_U + Rmom_Ust*Ust_U;
  Rmom_x = Rmom_Ust*Ust_x + 0.5*dF/(K^2*th0^2)*K_x;
  Rmom_CKx = m0;
  
  % shape equation
  h0 = 2*Hss/Hs+1-H; % at stagnation
  h0_Ust = 2*Hss_Ust/Hs - 2*Hss/Hs^2*Hs_Ust - H_Ust;
  [Hs2, Hs2_U2] = get_Hs(U(:,2), param);
  [Hs3, Hs3_U3] = get_Hs(U(:,3), param);
  dHs = Hs3-Hs2; dHs_U = [Hs2_U2, Hs3_U3];
  [Hss2, Hss2_U2] = get_Hss(U(:,2), param);
  [Hss3, Hss3_U3] = get_Hss(U(:,3), param);
  h2 = 2*Hss2/Hs2+1-H2;
  h2_U2 = 2*Hss2_U2/Hs2 - 2*Hss2/Hs2^2*Hs2_U2 - H2_U2;
  h3 = 2*Hss3/Hs3+1-H3;
  h3_U3 = 2*Hss3_U3/Hs3 - 2*Hss3/Hs3^2*Hs3_U3 - H3_U3;
  dh = h3-h2; dh_U = [-h2_U2, h3_U3];
  [D2, D2_U2] = get_cdutstag(U(:,2), param);
  [D3, D3_U3] = get_cdutstag(U(:,3), param);
  G2 = D2-0.5*F2; G2_U2 = D2_U2-0.5*F2_U2;
  G3 = D3-0.5*F3; G3_U3 = D3_U3-0.5*F3_U3;
  dG = G3-G2; dG_U = [-G2_U2, G3_U3];
  Rshape = th0/Hs*dHs + dh*th0 + h0*(CKx*th0 + 2*dth) - dG/(K*th0);
  Rshape_Ust = (th0_Ust-th0/Hs*Hs_Ust)/Hs*dHs + dh*th0_Ust + ...
      h0_Ust*(CKx*th0+2*dth) + h0*(CKx*th0_Ust) + dG/(K*th0^2)*th0_Ust;
  Rshape_U = th0/Hs*dHs_U + dh_U*th0 + h0*2*dth_U ...
      - dG_U/(K*th0) + dG/(K^2*th0)*K_U + Rshape_Ust*Ust_U;
  Rshape_x = Rshape_Ust*Ust_x + dG/(K^2*th0)*K_x;
  Rshape_CKx = h0*th0;
  
  % amplification (none, difference this time)
  Ramp = U(3,3)-U(3,2);
  Ramp_U = [0,0,-1,0, 0,0,1,0];
  Ramp_x = [0,0];
  Ramp_CKx = 0;
  
  % put into the global system
  Ig = 3*Ist(2) + (-2:0); Jg = [4*I(2)+(-3:0), 4*I(3)+(-3:0)];
  M.glob.R(Ig) = [Rmom; Rshape; Ramp];
  M.glob.R_U(Ig,:) = 0; M.glob.R_x(Ig,:) = 0; % clear out rows
  M.glob.R_U(Ig,Jg) = [Rmom_U; Rshape_U; Ramp_U];
  M.glob.R_x(Ig,Ist) = [Rmom_x; Rshape_x; Ramp_x];
  
  Jgue = 4*I;
  R_CKx = [Rmom_CKx; Rshape_CKx; Ramp_CKx];
  M.glob.R_U(Ig,Jgue) = M.glob.R_U(Ig,Jgue) + R_CKx*CKx_ue;
  M.glob.R_x(Ig,I) = M.glob.R_x(Ig,I) + R_CKx*CKx_x;
  
end


%-------------------------------------------------------------------------------
function system_test(M)
  
  M.oper.alpha = 3; % angle, in degrees
  M.oper.Ma = 0.0; % Mach number
  M.oper.Re = 5e5; % Reynolds number
  M.oper.viscous = true; % tests are viscous
  init_thermo(M); % thermodynamics
  
  N = 20;
  %Is = M.vsol.Is{1}; Is = Is(1:N); xi = M.isol.xi(Is); 
  %ds = M.glob.U(2,Is); ue = M.glob.U(4,Is);
  %xi = space_geom(8e-4, 0.042, N+1); xi = xi(2:end);
  xi = linspace(0,.042,N+1); xi = xi(2:end);
  th = ones(N,1)*1e-4; % initial theta
  ds = 2e-4 + 5e-4*xi - 1e-3*xi.^2; ds = ds'; % delta*
  
  ue = 1 - exp(-xi/.03); ue = ue';
  param = build_param(M, 1); % parameters

  nNewton = 100;
  for iNewton=1:nNewton
    
    % residual vector and Jacobian matrix
    R = zeros(2*N,1); R_U = sparse(2*N,2*N);
    
    % loop over stations
    for i=1:N
      im = i-1; if (im<1), im=1; end
      param.simi = (i==1);
      
      if (i==1)
        U1 = [th(i  ); ds(i  ); 0; ue(i  )];
        U2 = [th(i+1); ds(i+1); 0; ue(i+1)];
        x1 = xi(1); x2 = xi(2);
        dx = x2-x1;
        w1 =  x2/dx; w2 = -x1/dx;
        Ust = U1*w1 + U2*w2;
        
        %ff = 1.00;
        %K = ff*U1(4)/x1; % ue = K*x at stag
        %K_U = [0,0,0,ff/x1, 0,0,0,0];
        
        den = x1^2*x2 - x2^2*x1;
        wk1 = -x2^2/den; wk2 = x1^2/den;
        K = wk1*U1(4) + wk2*U2(4);
        K_U = [0,0,0,wk1, 0,0,0,wk2];
        
        xst = 1e-6*M.geom.chord; % xst needs to be small but nonzero
        Ust(4) = K*xst;
        Ust_U = [w1*eye(3,4), w2*eye(3,4); K_U*xst];
        param.turb = false; param.simi = true;  % similarity station flag  
        [R1, R1_Ust, ~] = residual_station(param, [xst,xst], [Ust, Ust], zeros(1,2));
        param.simi = false;
        R1_Ust = R1_Ust(:,1:4) + R1_Ust(:,5:8);
        
        R1_U = R1_Ust*Ust_U; 
        
        %dx_x = [-1,1];
        %w1_x = -w1/dx*dx_x + [ 0,1]/dx;
        %w2_x = -w2/dx*dx_x + [-1,0]/dx;
        %K_x = [-ff*U1(4)/x1^2, 0];
        %Ust_x = [U1(1:3)*w1_x + U2(1:3)*w2_x; K_x*xst];
        %R1_x = R1_Ust*Ust_x;

        I = [2*i-1,2*i]; J = [2*i-1,2*i,2*i+1,2*i+2]; JJ = [1,2,5,6];
        R(I) = R1(1:2);
        R_U(I,J) = R_U(I,J) + R1_U(1:2,JJ);
                        
      else
        x = xi([im,i]);
        U1 = [th(im); ds(im); 0; ue(im)];
        U2 = [th(i ); ds(i ); 0; ue(i )];
        [Ri, Ri_U, ~] = residual_station(param, x, [U1,U2], zeros(1,2));
        I = [2*i-1,2*i]; J = [2*im-1,2*im,2*i-1,2*i]; JJ = [1,2,5,6];
        R(I) = Ri(1:2);
        R_U(I,J) = R_U(I,J) + Ri_U(1:2,JJ);
      end
      
    end
    
    condest(R_U)    
    fprintf(1, 'iNewton = %d: norm(R) = %.5e\n', iNewton, norm(R));
    
    % solve/update
    dtd = -R_U\R;
    dth = dtd(1:2:end);
    dds = dtd(2:2:end);
    omega = 0.1; %0.1+0.9*(1-exp(-iNewton/2));
    if (iNewton > 20), omega = 1.; end
    th = th + omega*dth;
    ds = ds + omega*dds;
  end
  
  figure(1); clf; plot(xi,th,'bo-'); ylabel('theta');
  figure(2); clf; plot(xi,ds,'ro-'); ylabel('ds');
  figure(3); clf; plot(xi,ds./th,'ks-'); ylabel('H');
  figure(4); clf; plot(xi,ue,'m^-'); ylabel('ue');
  
end
  
  

%-------------------------------------------------------------------------------
function [upw, upw_U] = get_upw(U1, U2, param)
% calculates a local upwind factor (0.5 = trap; 1 = BE) based on two states
% INPUT
%   U1,U2 : first/upwind and second/downwind states (4x1 each)
%   param : parameter structure
% OUTPUT
%   upw   : scalar upwind factor
%   upw_U : 1x8 linearization vector, [upw_U1, upw_U2]
% DETAILS
%   Used to ensure a stable viscous discretization
%   Decision to upwind is made based on the shape factor change
  
  [Hk1, Hk1_U1] = get_Hk(U1, param);
  [Hk2, Hk2_U2] = get_Hk(U2, param);
  Z = zeros(size(Hk1_U1));
  Hut = 1.0; % triggering constant for upwinding
  C = 5.0; if (param.wake), C = 1.0; end
  Huc = C*Hut/Hk2^2; % only depends on U2
  Huc_U = [Z, -2*Huc/Hk2*Hk2_U2];
  aa = (Hk2-1)/(Hk1-1); sga = sign(aa);
  la = log(sga*aa);
  la_U = [-1/(Hk1-1)*Hk1_U1, 1/(Hk2-1)*Hk2_U2];
  Hls = la^2; Hls_U = 2*la*la_U;
  if (Hls > 15), Hls = 15; Hls_U = Hls_U*0; end
  upw = 1 - 0.5*exp(-Hls*Huc);
  upw_U = -0.5*exp(-Hls*Huc)*(-Hls_U*Huc-Hls*Huc_U);

end
  

%-------------------------------------------------------------------------------
function [f, f_U] = upwind(upw, upw_U, f1, f1_U1, f2, f2_U2)
% calculates an upwind average (and derivatives) of two scalars
% INPUT
%   upw, upw_U : upwind scalar and its linearization w.r.t. U1,U2
%   f1, f1_U   : first scalar and its linearization w.r.t. U1
%   f2, f2_U   : second scalar and its linearization w.r.t. U2
% OUTPUT
%   f    : averaged scalar
%   f_U  : linearization of f w.r.t. both states, [f_U1, f_U2]
  
  f = (1-upw)*f1 + upw*f2;
  f_U = (-upw_U)*f1 + upw_U*f2 + [(1-upw)*f1_U1, upw*f2_U2];
  
end


%-------------------------------------------------------------------------------
function [Hkc, rd] = slimit_Hkc(Hkc0)
% smooth limit of Hkc = function of Hk and Ret
% INPUT
%   Hkc0 : baseline value of Hkc
% OUTPUT
%   Hkc : smoothly limited value in defined range
%   rd  : derivative of Hkc w.r.t. the input Hkc0
  
% TODO: ping me
  Hl = .01; Hh = .05;
  if (Hkc0 < Hh)
    rn = (Hkc0-Hl)/(Hh-Hl); rn_Hkc0 = 1/(Hh-Hl);
    if (rn<0), rn=0.; rn_Hkc0 = 0.; end
    rf = 3*rn^2 - 2*rn^3; rf_rn = 6*rn - 6*rn^2;
    Hkc = Hl + rf*(Hh-Hl); rd = rf_rn*rn_Hkc0*(Hh-Hl);
  else
    Hkc = Hkc0; rd = 1;
  end
  
end


%-------------------------------------------------------------------------------
function [uq, uq_U] = get_uq(ds, ds_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param)
% calculates the equilibrium 1/ue*due/dx
% INPUT
%   ds, ds_U   : delta star and linearization (1x4)
%   cf, cf_U   : skin friction and linearization (1x4)
%   Hk, Hk_U   : kinematic shape parameter and linearization (1x4)
%   Ret, Ret_U : theta Reynolds number and linearization (1x4)
%   param      : parameter structure
% OUTPUT
%   uq, uq_U   : equilibrium 1/ue*due/dx and linearization w.r.t. state (1x4)

  beta = param.GB; A = param.GA; C = param.GC;
  if (param.wake), A = A*param.Dlr; C = 0; end
  % limit Hk (TODO smooth/eliminate)
  if (param.wake) && (Hk < 1.00005), Hk = 1.00005; Hk_U = Hk_U*0; end
  if (~param.wake) && (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
  Hkc = Hk - 1 - C/Ret;
  Hkc_U = Hk_U + C/Ret^2*Ret_U;
  %[Hkc, rd] = slimit_Hkc(Hkc); Hkc_U = rd*Hkc_U; % smooth limiting of Hkc  
  
  if (Hkc < .01), Hkc = .01; Hkc_U = Hkc_U*0.; end
  ut = 0.5*cf - (Hkc/(A*Hk))^2;
  ut_U = 0.5*cf_U - 2*(Hkc/(A*Hk))*(Hkc_U/(A*Hk) - Hkc/(A*Hk^2)*Hk_U);
  uq = ut/(beta*ds);
  uq_U = ut_U/(beta*ds) - uq/ds * ds_U;
  
end


%-------------------------------------------------------------------------------
function [cttr, cttr_U] = get_cttr(U, param)
% calculates root of the shear stress coefficient at transition
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cttr, cttr_U : sqrt(shear stress coeff) and its lin w.r.t. U (1x4)
% DETAILS
%   used to initialize the first turb station after transition
  
  param.wake = false; % transition happens just before the wake starts
  [cteq, cteq_U] = get_cteq(U, param);
  [Hk, Hk_U] = get_Hk(U, param);
  if (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
  C = param.CtauC; E = param.CtauE;
  c = C*exp(-E/(Hk-1)); c_U = c*E/(Hk-1)^2*Hk_U;
  cttr = c*cteq; cttr_U = c_U*cteq + c*cteq_U;

end


%-------------------------------------------------------------------------------
function [cteq, cteq_U] = get_cteq(U, param)
% calculates root of the equilibrium shear stress coefficient: sqrt(ctau_eq)
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cteq, cteq_U : sqrt(equilibrium shear stress) and its lin w.r.t. U (1x4)
% DETAILS
%   uses equilibrium shear stress correlations
  CC = 0.5/(param.GA^2*param.GB); C = param.GC;
  [Hk, Hk_U] = get_Hk(U, param);
  [Hs, Hs_U] = get_Hs(U, param);
  [H, H_U] = get_H(U);
  [Ret, Ret_U] = get_Ret(U, param);
  [Us, Us_U] = get_Us(U, param);
  if (param.wake)
    if (Hk < 1.00005), Hk = 1.00005; Hk_U = Hk_U*0; end
    Hkc = Hk - 1; Hkc_U = Hk_U;
  else
    if (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
    Hkc = Hk - 1 - C/Ret;
    Hkc_U = Hk_U + C/Ret^2*Ret_U;
    %[Hkc, rd] = slimit_Hkc(Hkc); Hkc_U = rd*Hkc_U; % smooth limiting of Hkc  
    if (Hkc < 0.01), Hkc = 0.01; Hkc_U = Hkc_U*0.; end
  end
  num = CC*Hs*(Hk-1)*Hkc^2;
  num_U = CC*(Hs_U*(Hk-1)*Hkc^2 + Hs*Hk_U*Hkc^2 + Hs*(Hk-1)*2*Hkc*Hkc_U);
  den = (1-Us)*H*Hk^2;
  den_U = (-Us_U)*H*Hk^2 + (1-Us)*H_U*Hk^2 + (1-Us)*H*2*Hk*Hk_U;
  cteq = sqrt(num/den);
  cteq_U = 0.5/cteq*(num_U/den - num/den^2*den_U);

end


%-------------------------------------------------------------------------------
function [Hs, Hs_U] = get_Hs(U, param)
% calculates Hs = Hstar = K.E. shape parameter, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   Hs, Hs_U : Hstar and its lin w.r.t. U (1x4)
% DETAILS
%   Hstar is the ratio theta*/theta, where theta* is the KE thicknes
  [Hk, Hk_U] = get_Hk(U, param);
  
  % limit Hk (TODO smooth/eliminate)
  if (param.wake) && (Hk < 1.00005), Hk = 1.00005; Hk_U = Hk_U*0; end
  if (~param.wake) && (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
  
  if (param.turb) % turbulent
    Hsmin = 1.5; dHsinf = .015;
    [Ret, Ret_U] = get_Ret(U, param);
    % limit Re_theta and dependence
    Ho = 4; Ho_U = 0.;
    if (Ret > 400), Ho = 3 + 400/Ret; Ho_U = -400/Ret^2*Ret_U; end
    Reb = Ret; Reb_U = Ret_U;
    if (Ret < 200), Reb = 200; Reb_U = Reb_U*0; end
    if (Hk < Ho)  % attached branch
      Hr = (Ho-Hk)/(Ho-1);
      Hr_U = (Ho_U - Hk_U)/(Ho-1) - (Ho-Hk)/(Ho-1)^2*Ho_U;
      aa = (2-Hsmin-4/Reb)*Hr^2;
      aa_U = (4/Reb^2*Reb_U)*Hr^2 + (2-Hsmin-4/Reb)*2*Hr*Hr_U;
      Hs = Hsmin + 4/Reb + aa * 1.5/(Hk+.5);
      Hs_U = -4/Reb^2*Reb_U + aa_U*1.5/(Hk+.5) - aa*1.5/(Hk+.5)^2*Hk_U;
    else % separated branch
      lrb = log(Reb); lrb_U = 1/Reb*Reb_U;
      aa = Hk - Ho + 4/lrb;
      aa_U = Hk_U - Ho_U - 4/lrb^2*lrb_U;
      bb = .007*lrb/aa^2 + dHsinf/Hk;
      bb_U = .007*(lrb_U/aa^2 - 2*lrb/aa^3*aa_U) - dHsinf/Hk^2*Hk_U;
      Hs = Hsmin + 4/Reb + (Hk-Ho)^2*bb;
      Hs_U = -4/Reb^2*Reb_U + 2*(Hk-Ho)*(Hk_U-Ho_U)*bb + (Hk-Ho)^2*bb_U;
    end
    % slight Mach number correction
    [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
    den = 1+.014*M2; den_M2 = .014;
    Hs = (Hs+.028*M2)/den;
    Hs_U = (Hs_U+.028*M2_U)/den - Hs/den*den_M2*M2_U;
  else % laminar
    a = Hk-4.35;
    if (Hk < 4.35)
      num = .0111*a^2 - .0278*a^3;
      Hs = num/(Hk+1) + 1.528 - .0002*(a*Hk)^2;
      Hs_Hk = (.0111*2*a - .0278*3*a^2)/(Hk+1) - num/(Hk+1)^2 - .0002*2*a*Hk*(Hk+a);
    else
      Hs = .015*a^2/Hk + 1.528;
      Hs_Hk = .015*2*a/Hk - .015*a^2/Hk^2;
    end
    Hs_U = Hs_Hk*Hk_U;
  end

end


%-------------------------------------------------------------------------------
function [cp, cp_u] = get_cp(u, param)
% calculates pressure coefficient from speed, with compressibility correction
% INPUT
%   u     : speed
%   param : parameter structure
% OUTPUT
%   cp, cp_U : pressure coefficient and its linearization w.r.t. u
% DETAILS
%   Karman-Tsien correction is included

  Vinf = param.Vinf;
  cp = 1-(u/Vinf).^2; cp_u = -2*u/Vinf^2;
  if (param.Minf > 0)
    l = param.KTl; b = param.KTb;
    den = b+0.5*l*(1+b)*cp; den_cp = 0.5*l*(1+b);
    cp = cp./den; cp_u = cp_u .* (1-cp*den_cp)./den;
  end
  
end


%-------------------------------------------------------------------------------
function [uk, uk_u] = get_uk(u, param)
% calculates Karman-Tsien corrected speed
% INPUT
%   u     : incompressible speed
%   param : parameter structure
% OUTPUT
%   uk, uk_u : compressible speed and its linearization w.r.t. u
% DETAILS
%   Uses the Karman-Tsien correction, Minf from param
  
  if (param.Minf > 0)
    l = param.KTl; Vinf = param.Vinf;
    den = 1-l*(u/Vinf).^2; den_u = -2*l*u/Vinf^2;
    uk = u*(1-l)./den; uk_u = (1-l)./den - (uk./den).*den_u;
  else 
    uk = u; uk_u = 1; 
  end
  
end


%-------------------------------------------------------------------------------
function [M2, M2_U] = get_Mach2(U, param)
% calculates squared Mach number
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   M2, M2_U : squared Mach number and its linearization w.r.t. U (1x4)
% DETAILS
%   Uses constant total enthalpy from param.H0
%   The speed of sound varies; depends on enthalpy, which depends on speed
%   The compressible edge speed must be used
  
  if (param.Minf > 0)
    H0 = param.H0; g = param.gam;
    [uk, uk_u] = get_uk(U(4), param);
    c2 = (g-1)*(H0-0.5*uk^2); c2_uk = (g-1)*(-uk); % squared speed of sound
    M2 = uk^2/c2; M2_uk = 2*uk/c2 - M2/c2*c2_uk; M2_U = [0,0,0,M2_uk*uk_u];
  else
    M2 = 0.; M2_U = zeros(1,4); 
  end
  
end


%-------------------------------------------------------------------------------
function [H, H_U] = get_H(U)
% calculates H = shape parameter = delta*/theta, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   H, H_U : shape parameter and its linearization w.r.t. U (1x4)
% DETAILS
%   H is the ratio of the displacement thickness to the momentum thickness
%   In U, the ds entry should be (delta*-wgap) ... i.e wake gap taken out
%   When the real H is needed with wake gap, Hw is calculated and added
  
  H = U(2)/U(1);
  H_U = [-H/U(1), 1/U(1), 0, 0];
  
end


%-------------------------------------------------------------------------------
function [Hw, Hw_U] = get_Hw(U, wgap)
% calculates Hw = wake gap shape parameter = wgap/theta
% INPUT
%   U    : state vector [th; ds; sa; ue]
%   wgap : wake gap
% OUTPUT
%   Hw, Hw_U : wake gap shape parameter and its linearization w.r.t. U (1x4)
% DETAILS
%   Hw is the ratio of the wake gap to the momentum thickness
%   The wake gap is the TE gap extrapolated into the wake (dead air region)
  
  Hw = wgap/U(1); % wgap/th
  Hw_U = [-Hw/U(1),0,0,0];
  
end


%-------------------------------------------------------------------------------
function [Hk, Hk_U] = get_Hk(U, param)
% calculates Hk = kinematic shape parameter, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   Hk, Hk_U : kinematic shape parameter and its linearization w.r.t. U (1x4)
% DETAILS
%   Hk is like H but with no density in the integrals defining th and ds
%   So it is exactly the same when density is constant (= freestream)
%   Here, it is computed from H with a correlation using the Mach number
  
  [H, H_U] = get_H(U);
  
  if (param.Minf > 0)
    [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
    den = (1+0.113*M2); den_M2 = 0.113;
    Hk = (H-0.29*M2)/den;
    Hk_U = (H_U-0.29*M2_U)/den - Hk/den*den_M2*M2_U;
  else
    Hk = H; Hk_U = H_U;
  end
  
end


%-------------------------------------------------------------------------------
function [Hss, Hss_U] = get_Hss(U, param)
% calculates Hss = density shape parameter, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   Hss, Hss_U : density shape parameter and its linearization w.r.t. U (1x4)
% DETAILS
  
  [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
  [Hk, Hk_U] = get_Hk(U,param);
  num = 0.064/(Hk-0.8) + 0.251; num_U = -.064/(Hk-0.8)^2*Hk_U;
  Hss = M2*num; Hss_U = M2_U*num + M2*num_U;
  
end


%-------------------------------------------------------------------------------
function [de, de_U] = get_de(U, param)
% calculates simplified BL thickness measure
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   de, de_U : BL thickness "delta" and its linearization w.r.t. U (1x4)
% DETAILS
%   delta is delta* incremented with a weighted momentum thickness, theta
%   The weight on theta depends on Hk, and there is an overall cap
  
  [Hk, Hk_U] = get_Hk(U, param);
  aa = 3.15 + 1.72/(Hk-1); aa_U = -1.72/(Hk-1)^2*Hk_U;
  de = U(1)*aa + U(2); de_U = [aa,1,0,0] + U(1)*aa_U;
  dmx = 12.0;
  if (de > dmx*U(1)), de = dmx*U(1); de_U = [dmx,0,0,0]; end
  
end


%-------------------------------------------------------------------------------
function [rho, rho_U] = get_rho(U, param)
% calculates the density (useful if compressible)
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   rho, rho_U : density and linearization
% DETAILS
%   If compressible, rho is calculated from stag rho + isentropic relations
  
  if (param.Minf > 0)
    [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
    [uk, uk_u] = get_uk(U(4), param); % corrected speed
    H0 = param.H0; gmi = param.gam-1;
    den = 1+0.5*gmi*M2; den_M2 = 0.5*gmi;
    rho = param.rho0/den^(1/gmi); rho_U = (-1/gmi)*rho/den*den_M2*M2_U;
  else
    rho = param.rho0; 
    rho_U = zeros(1,4);
  end
  
end


%-------------------------------------------------------------------------------
function [Ret, Ret_U] = get_Ret(U, param)
% calculates theta Reynolds number, Re_theta, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   Ret, Ret_U : Reynolds number based on the momentum thickness, linearization
% DETAILS
%   Re_theta = rho*ue*theta/mu
%   If compressible, rho is calculated from stag rho + isentropic relations
%   ue is the edge speed and must be comressibility corrected
%   mu is the dynamic viscosity, from Sutherland's law if compressible
  
  if (param.Minf > 0)
    [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
    [uk, uk_u] = get_uk(U(4), param); % corrected speed
    H0 = param.H0; gmi = param.gam-1; Ts = param.Tsrat;
    Tr = 1-0.5*uk^2/H0; Tr_uk = -uk/H0; % edge/stagnation temperature ratio
    f = Tr^1.5*(1+Ts)/(Tr+Ts); f_Tr = 1.5*f/Tr-f/(Tr+Ts); % Sutherland's ratio
    mu = param.mu0*f; mu_uk = param.mu0*f_Tr*Tr_uk; % local dynamic viscosity
    den = 1+0.5*gmi*M2; den_M2 = 0.5*gmi;
    rho = param.rho0/den^(1/gmi); rho_U = (-1/gmi)*rho/den*den_M2*M2_U; % density
    Ret = rho*uk*U(1)/mu;
    Ret_U = rho_U*uk*U(1)/mu + (rho*U(1)/mu-Ret/mu*mu_uk)*[0,0,0,uk_u] + rho*uk/mu*[1,0,0,0];
  else
    Ret = param.rho0*U(1)*U(4)/param.mu0;
    Ret_U = [U(4), 0, 0, U(1)]/param.mu0;
  end
  
end


%-------------------------------------------------------------------------------
function [cf, cf_U] = get_cf(U, param)
% calculates cf = skin friction coefficient, from U
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cf, cf_U : skin friction coefficient and its linearization w.r.t. U (1x4)
% DETAILS
%   cf is the local skin friction coefficient = tau/(0.5*rho*ue^2)
%   Correlations are used based on Hk and Re_theta

  if (param.wake), cf = 0; cf_U = zeros(1,4); return; end % zero cf in wake
  [Hk, Hk_U] = get_Hk(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  
  % TODO: limit Hk
  
  if (param.turb) % turbulent cf
    [M2, M2_U] = get_Mach2(U, param); % squared edge Mach number
    Fc = sqrt(1+0.5*(param.gam-1)*M2);
    Fc_U = 0.5/Fc*0.5*(param.gam-1)*M2_U;
    aa = -1.33*Hk; aa_U = -1.33*Hk_U;
    %if (aa < -20), aa = -20; aa_U = aa_U*0; warning('aa in cfturb'); end
    % smooth limiting of aa
    if (aa < -17), aa = -20+3*exp((aa+17)/3); aa_U = (aa+20)/3*aa_U; end % TODO: ping me  
    bb = log(Ret/Fc); bb_U = Ret_U/Ret - Fc_U/Fc;
    if (bb < 3), bb = 3; bb_U = bb_U*0; end
    bb = bb/log(10); bb_U = bb_U/log(10);
    cc = -1.74 - 0.31*Hk; cc_U = -0.31*Hk_U;
    dd = tanh(4.0-Hk/0.875); dd_U = (1-dd^2)*(-Hk_U/0.875);
    cf0 = 0.3*exp(aa)*bb^cc;
    cf0_U = cf0*aa_U + 0.3*exp(aa)*cc*bb^(cc-1)*bb_U + cf0*log(bb)*cc_U;
    cf = (cf0 + 1.1e-4*(dd-1))/Fc;
    cf_U = (cf0_U + 1.1e-4*dd_U)/Fc - cf/Fc*Fc_U;
  else % laminar cf
    if (Hk < 5.5)
      num = .0727*(5.5-Hk)^3/(Hk+1) - .07;
      num_Hk = .0727*(3*(5.5-Hk)^2/(Hk+1)*(-1) - (5.5-Hk)^3/(Hk+1)^2);
    else
      num = .015*(1-1./(Hk-4.5))^2 - .07;
      num_Hk = .015*2*(1-1./(Hk-4.5))/(Hk-4.5)^2;
    end
    cf = num/Ret;
    cf_U = num_Hk/Ret*Hk_U - num/Ret^2*Ret_U;
  end

end


%-------------------------------------------------------------------------------
function [cfxt, cfxt_U, cfxt_x] = get_cfxt(U, x, param)
% calculates cf*x/theta from the state
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   x     : distance along wall (xi)
%   param : parameter structure
% OUTPUT
%   cfxt,  : the combination cf*x/theta (calls cf function)
%   cfxt_U : linearization w.r.t. U (1x4)
%   cfxt_x : linearization w.r.t x (scalar)  
% DETAILS
%   This combination appears in the momentum and shape parameter equations
  
  [cf, cf_U] = get_cf(U, param);
  cfxt = cf*x/U(1); 
  cfxt_U = cf_U*x/U(1); cfxt_U(1) = cfxt_U(1) - cfxt/U(1);
  cfxt_x = cf/U(1);
  
end


%-------------------------------------------------------------------------------
function [F, F_U] = get_cfutstag(U, param)
% calculates cf*ue*theta, used in stagnation station calculations
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   F, F_U : value and linearization of cf*ue*theta
% DETAILS
%   Only for stagnation and laminar

  U(4) = 0; [Hk, Hk_U] = get_Hk(U, param);

  if (Hk < 5.5)
    num = .0727*(5.5-Hk)^3/(Hk+1) - .07;
    num_Hk = .0727*(3*(5.5-Hk)^2/(Hk+1)*(-1) - (5.5-Hk)^3/(Hk+1)^2);
  else
    num = .015*(1-1./(Hk-4.5))^2 - .07;
    num_Hk = .015*2*(1-1./(Hk-4.5))/(Hk-4.5)^2;
  end
  nu = param.mu0/param.rho0;
  F = nu*num;
  F_U = nu*num_Hk*Hk_U;

end


%-------------------------------------------------------------------------------
function [D, D_U] = get_cdutstag(U, param)
% calculates cDi*ue*theta, used in stagnation station calculations
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   D, D_U : value and linearization of cDi*ue*theta
% DETAILS
%   Only for stagnation and laminar

  U(4) = 0; [Hk, Hk_U] = get_Hk(U, param);
  
  if (Hk<4)
    num = .00205*(4-Hk)^5.5 + .207;
    num_Hk = .00205*5.5*(4-Hk)^4.5*(-1);
  else
    Hk1 = Hk-4;
    num = -.0016*Hk1^2/(1+.02*Hk1^2) + .207;
    num_Hk = -.0016*(2*Hk1/(1+.02*Hk1^2) - Hk1^2/(1+.02*Hk1^2)^2*.02*2*Hk1);
  end
  
  nu = param.mu0/param.rho0;
  D = nu*num;
  D_U = nu*num_Hk*Hk_U;

end


%-------------------------------------------------------------------------------
function [cDixt, cDixt_U, cDixt_x] = get_cDixt(U, x, param)
% calculates cDi*x/theta from the state
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   x     : distance along wall (xi)
%   param : parameter structure
% OUTPUT
%   cDixt,  : the combination cDi*x/theta (calls cDi function)
%   cDixt_U : linearization w.r.t. U (1x4)
%   cDixt_x : linearization w.r.t x (scalar)  
% DETAILS
%   cDi is the dissipation function
  
  [cDi, cDi_U] = get_cDi(U, param);
  cDixt = cDi*x/U(1); 
  cDixt_U = cDi_U*x/U(1); cDixt_U(1) = cDixt_U(1) - cDixt/U(1);
  cDixt_x = cDi/U(1);
  
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi(U, param)
% calculates cDi = dissipation function = 2*cD/H*, from the state
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   cD is the dissipation coefficient, int(tau*du/dn*dn)/(rho*ue^3)
%   The combination with H* appears in the shape parameter equation
  
  if (param.turb) % turbulent includes wake
    
    % initialize to 0; will add components that are needed
    cDi = 0; cDi_U = zeros(1,4);
    
    if (~param.wake)
      % turbulent wall contribution (0 in the wake) 
      [cDi0, cDi0_U] = get_cDi_turbwall(U, param);
      cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U;
      [cDil, cDil_U] = get_cDi_lam(U, param); % for max check
    else
      [cDil, cDil_U] = get_cDi_lamwake(U, param); % for max check
    end
    
    % outer layer contribution
    [cDi0, cDi0_U] = get_cDi_outer(U, param);
    cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U;
    
    % laminar stress contribution
    [cDi0, cDi0_U] = get_cDi_lamstress(U, param);
    cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U;
        
    % maximum check
    if (cDil > cDi), cDi = cDil; cDi_U = cDil_U; end
    
    % double dissipation in the wake
    if (param.wake), cDi = 2*cDi; cDi_U = 2*cDi_U; end
  else
    % just laminar dissipation
    [cDi, cDi_U] = get_cDi_lam(U, param);
  end
  
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi_turbwall(U, param)
% calculates the turbulent wall contribution to cDi
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   This is one contribution to the dissipation function cDi = 2*cD/H*
  
  cDi = 0; cDi_U = zeros(1,4); if (param.wake), return; end
  
  % get cf, Hk, Hs, Us
  [cf, cf_U] = get_cf(U, param);
  [Hk, Hk_U] = get_Hk(U, param);
  [Hs, Hs_U] = get_Hs(U, param);
  [Us, Us_U] = get_Us(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  
  lr = log(Ret); lr_U = Ret_U/Ret;
  Hmin = 1 + 2.1/lr; Hmin_U = -2.1/lr^2*lr_U;
  aa = tanh((Hk-1)/(Hmin-1)); fac = 0.5 + 0.5*aa;
  fac_U = 0.5*(1-aa^2)*(Hk_U/(Hmin-1)-(Hk-1)/(Hmin-1)^2*Hmin_U);

  cDi = 0.5*cf*Us*(2/Hs)*fac;
  cDi_U = cf_U*Us/Hs*fac + cf*Us_U/Hs*fac - cDi/Hs*Hs_U + cf*Us/Hs*fac_U;
    
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi_lam(U, param)
% calculates the laminar dissipation function cDi
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   This is one contribution to the dissipation function cDi = 2*cD/H*
  
  % first get Hk and Ret
  [Hk, Hk_U] = get_Hk(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  
  if (Hk<4)
    num = .00205*(4-Hk)^5.5 + .207;
    num_Hk = .00205*5.5*(4-Hk)^4.5*(-1);
  else
    Hk1 = Hk-4;
    num = -.0016*Hk1^2/(1+.02*Hk1^2) + .207;
    num_Hk = -.0016*(2*Hk1/(1+.02*Hk1^2) - Hk1^2/(1+.02*Hk1^2)^2*.02*2*Hk1);
  end
  cDi = num/Ret;
  cDi_U = num_Hk/Ret*Hk_U - num/Ret^2*Ret_U;
  
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi_lamwake(U, param)
% laminar wake dissipation function cDi
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   This is one contribution to the dissipation function cDi = 2*cD/H*
  
  param.turb = false; % force laminar
  
  % dependencies
  [Hk, Hk_U] = get_Hk(U, param);
  [Hs, Hs_U] = get_Hs(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  HsRet = Hs*Ret;
  HsRet_U = Hs_U*Ret + Hs*Ret_U;
  
  num = 2*1.1*(1-1/Hk)^2*(1/Hk);
  num_Hk = 2*1.1*(2*(1-1/Hk)*(1/Hk^2)*(1/Hk)+(1-1/Hk)^2*(-1/Hk^2));
  cDi = num/HsRet;
  cDi_U = num_Hk*Hk_U/HsRet - num/HsRet^2*HsRet_U;
      
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi_outer(U, param)
% turbulent outer layer contribution to dissipation function cDi
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   This is one contribution to the dissipation function cDi = 2*cD/H*
  
  if (~param.turb), cDi=0; cDi_U = zeros(1,4); return; end % for pinging
  
  % first get Hs, Us
  [Hs, Hs_U] = get_Hs(U, param);
  [Us, Us_U] = get_Us(U, param);

  % shear stress: note, state stores ct^.5
  ct = U(3)^2; ct_U = [0,0,2*U(3),0];
  
  cDi = ct*(0.995-Us)*2/Hs;
  cDi_U = ct_U*(0.995-Us)*2/Hs + ct*(-Us_U)*2/Hs - ct*(0.995-Us)*2/Hs^2*Hs_U;
  
end


%-------------------------------------------------------------------------------
function [cDi, cDi_U] = get_cDi_lamstress(U, param)
% laminar stress contribution to dissipation function cDi
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
% DETAILS
%   This is one contribution to the dissipation function cDi = 2*cD/H*
  
  % first get Hs, Us, and Ret
  [Hs, Hs_U] = get_Hs(U, param);
  [Us, Us_U] = get_Us(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  HsRet = Hs*Ret;
  HsRet_U = Hs_U*Ret + Hs*Ret_U;
  
  num = 0.15*(0.995-Us)^2*2;
  num_Us = 0.15*2*(0.995-Us)*(-1)*2;
  cDi = num/HsRet;
  cDi_U = num_Us*Us_U/HsRet - num/HsRet^2*HsRet_U;
  
end


%-------------------------------------------------------------------------------
function [Us, Us_U] = get_Us(U, param)
% calculates the normalized wall slip velocity Us
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   Us, Us_U : normalized wall slip velocity and its linearization w.r.t. U (1x4)
  
  [Hs, Hs_U] = get_Hs(U, param);
  [Hk, Hk_U] = get_Hk(U, param);
  [H, H_U] = get_H(U);
  
  % limit Hk (TODO smooth/eliminate)
  if (param.wake) && (Hk < 1.00005), Hk = 1.00005; Hk_U = Hk_U*0; end
  if (~param.wake) && (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
  
  beta = param.GB; bi = 1./beta;
  Us = 0.5*Hs*(1-bi*(Hk-1)/H);
  Us_U = 0.5*Hs_U*(1-bi*(Hk-1)/H) + 0.5*Hs*(-bi*(Hk_U)/H +bi*(Hk-1)/H^2*H_U);  
  % limits
  if (~param.wake && (Us>0.95   )), Us = 0.98   ; Us_U = Us_U*0; end
  if ( param.wake && (Us>0.99995)), Us = 0.99995; Us_U = Us_U*0; end  
  
end  


%-------------------------------------------------------------------------------
function [damp, damp_U] = get_damp(U, param)
% calculates the amplification rate, dn/dx, used in predicting transition
% INPUT
%   U     : state vector [th; ds; sa; ue]
%   param : parameter structure
% OUTPUT
%   damp, damp_U : amplification rate and its linearization w.r.t. U (1x4)
% DETAILS
%   damp = dn/dx is used in the amplification equation, prior to transition
  
  [Hk, Hk_U] = get_Hk(U, param);
  [Ret, Ret_U] = get_Ret(U, param);
  th = U(1);
    
  % limit Hk (TODO smooth/eliminate)
  if (Hk < 1.05), Hk = 1.05; Hk_U = Hk_U*0; end
  
  Hmi = 1/(Hk-1); Hmi_U = -Hmi^2*Hk_U;
  aa = 2.492*Hmi^0.43; aa_U = 0.43*aa/Hmi*Hmi_U;
  bb = tanh(14*Hmi-9.24); bb_U = (1-bb^2)*14*Hmi_U;
  lrc = aa + 0.7*(bb+1); lrc_U = aa_U + 0.7*bb_U;
  lten = log(10); lr = log(Ret)/lten; lr_U = (1/Ret)*Ret_U/lten;
  dl = .1;  % changed from .08 to make smoother
  damp = 0; damp_U = zeros(size(U')); % default no amplification
  if (lr >= lrc-dl)
    rn = (lr-(lrc-dl))/(2*dl); rn_U = (lr_U - lrc_U)/(2*dl);
    if (rn >= 1)
      rf = 1; rf_U = zeros(size(U'));
    else
      rf = 3*rn^2-2*rn^3; rf_U = (6*rn-6*rn^2)*rn_U;
    end
    ar = 3.87*Hmi-2.52; ar_U = 3.87*Hmi_U;
    ex = exp(-ar^2); ex_U = ex*(-2*ar*ar_U);
    da = 0.028*(Hk-1)-0.0345*ex; da_U = 0.028*Hk_U-0.0345*ex_U;
    af = -0.05+2.7*Hmi-5.5*Hmi^2+3*Hmi^3+0.1*exp(-20*Hmi);
    af_U = (2.7-11*Hmi+9*Hmi^2-1*exp(-20*Hmi))*Hmi_U;
    damp = rf*af*da/th;
    damp_U = (rf_U*af*da + rf*af_U*da + rf*af*da_U)/th - damp/th*[1,0,0,0];
  end
    
  % extra amplification to ensure dn/dx > 0 near ncrit
  ncrit = param.ncrit;
  
  Cea = 5; nx = Cea*(U(3)-ncrit); nx_U = Cea*[0,0,1,0];
  eex = 1+tanh(nx); eex_U = (1-tanh(nx)^2)*nx_U;
  
  ed = eex*.001/th;
  ed_U = eex_U*.001/th - ed/th*[1,0,0,0];
  damp = damp + ed;
  damp_U = damp_U + ed_U;
    
end  


%-------------------------------------------------------------------------------
function [E, rate] = check_ping(ep, v, v_u, sname)
% checks convergence of 3 values/derivatives
% INPUT
%   v     : cell array of three function evaluations at 0,+ep,+2*ep
%   v_u   : cell array of three derivative evaluations at 0,+ep,+2*ep
%   sname : descriptive name of where values came from for printing
% OUTPUT
%   E     : error values for two finite-difference comparisons
%   rate  : convergence rate, also printed
  
  E = zeros(1,2);
  for i = 1:2, E(i) = norm((v{1+i}-v{1})/(ep*i) - 0.5*(v_u{1} + v_u{1+i})); end
  rate = log2(E(2)/E(1));
  fprintf(1, '%s ping error convergence rate = %.4f\n', sname, rate);

end

  
%-------------------------------------------------------------------------------
function ping_test(M)
% checks derivatives of various functions through finite-difference pinging
% INPUT
%   M : mfoil class
% OUTPUT
%   printouts of rates (2 = second order expected).
  
  M.oper.alpha = 3; % angle, in degrees
  M.oper.Ma = 0.4; % Mach number
  M.oper.viscous = true; % tests are viscous
  rng(0); % for consistent pseudo random numbers
  M.param.verb = 2; % to minimize prints to screen
  
  % freestream Reynolds numbers
  Rev = [2e3, 1e5];
  
  % laminar/turbulent test states: th, ds, sa, ue
  Uv = [0.01, 0.02, 8.4, 0.9;  0.023, 0.05, .031, 1.1]';
  
  % functions to test
  fv = {@get_Hk, @get_Ret, @get_cf, @get_cDi, @get_Hs, @get_Us, ...
        @get_cDi_turbwall, @get_cDi_lam, @get_cDi_lamwake, @get_cDi_outer, ...
        @get_cDi_lamstress, @get_cteq, @get_cttr, @get_de, @get_damp, ...
        @get_Mach2, @get_Hss, @residual_station};
  
  % ping tests
  sturb = {'lam', 'turb', 'wake'};
  for iRe = 1:length(Rev) % loop over Reynolds numbers
    M.oper.Re = Rev(iRe); 
    init_thermo(M);
    param = build_param(M, 1); 
    for it = 1:3 % loop over lam, turb, wake
      param.turb = (it>1); param.wake = (it==3);
      for ih = 1:length(fv) % loop over functions
        U = Uv(:,min(it,2)); srates = ''; smark = ''; serr = ''; f = fv{ih};
        if (isequal(f,@residual_station)), U = [U; U.*[1.1,.8,.9,1.2]']; end
        for k = 1:length(U) % test all state component derivatives
          ep = 1e-4; E = zeros(1,2); 
          if (isequal(f,@residual_station))
            xi = [0.7,0.8]; Aux = [.002, .0018]; dx = [-.2,.3];
            [v0, v_U0, v_x0] = f(M,param,xi,reshape(U,4,2),Aux);
            for iep = 1:2 % test with two epsilons
              U(k) = U(k) + ep; xi = xi + ep*dx;
              [v1, v_U1, v_x1] = f(M,param,xi,reshape(U,4,2),Aux);
              U(k) = U(k) - ep; xi = xi - ep*dx;
              E(iep) = norm((v1-v0)/ep - 0.5*(v_U1(:,k) + v_U0(:,k) + (v_x0+v_x1)*dx'));
              ep = ep/2;
            end
          else
            [v0, v_U0] = f(U, param);
            for iep = 1:2 % test with two epsilons
              U(k) = U(k) + ep; [v1, v_U1] = f(U, param); U(k) = U(k) - ep;
              E(iep) = abs((v1-v0)/ep - 0.5*(v_U1(k) + v_U0(k)));
              ep = ep/2;
            end
          end
          srate = ' N/A'; skip = false; %isequal(f,@get_Ret);
          if (~skip) && (E(2)>5e-11)
            m = log2(E(1)/E(2)); srate = sprintf('%4.1f', m);
            if (m<1.5), smark = '<==='; end
          end
          srates = [srates, ' ', srate];
          serr = [serr, sprintf(' %.2e', E(2))];
        end  
        vprint(param, 0, '%-18s %-5s err=[%s]  rates=[%s] %s\n', ...
               func2str(f), sturb{it}, serr, srates, smark);
      end
    end
  end
  
  % transition residual ping
  M.oper.Re = 2e6; init_thermo(M); param = build_param(M,1);
  U = [0.01, 0.02, 8.95, 0.9;  0.013, 0.023, .028, 0.85]'; x = [0.7,0.8]'; Aux = [0,0]; 
  rng(0); dU = rand(size(U)); dx = rand(size(x)); ep = 1e-4; v = cell(1,3); v_u = v;
  for ie = 1:3
    [R, R_U, R_x] = residual_transition(M, param, x, U, Aux);
    v{ie} = R; v_u{ie} = R_U*reshape(dU,8,1) + R_x*dx;
    U = U + ep*dU; x = x + ep*dx;
  end
  check_ping(ep, v, v_u, 'transition residual');  

  % stagnation residual ping
  M.oper.Re = 1e6; M.oper.alpha = 1; init_thermo(M); param = build_param(M,1);
  U = [0.00004673616, 0.000104289, 0, 0.11977917547]'; x = 4.590816441485401e-05; Aux = [0,0]; 
  rng(0); dU = rand(size(U)); dx = rand(size(x)); ep = 1e-6; v = cell(1,3); v_u = v;
  for ie = 1:3
    param.simi = true;
    [R, R_U, R_x] = residual_station(param, [x,x], [U,U], Aux);
    param.simi = false;
    v{ie} = R; v_u{ie} = (R_U(:,1:4)+R_U(:,5:8))*dU + (R_x(:,1)+R_x(:,2))*dx;
    U = U + ep*dU; x = x + ep*dx;
  end
  check_ping(ep, v, v_u, 'stagnation residual');
  
  % stagnation state ping
  M.oper.Re = 5e5; init_thermo(M); param = build_param(M,1);
  U = [5e-5, 1.1e-4, 0, .03486; 4.9e-5, 1.09e-4, 0, .07397]'; x = [5.18e-4, 1.1e-3]';
  rng(0); dU = rand(size(U)); dx = rand(size(x)); ep = 1e-6; v = cell(1,3); v_u = v;
  for ie = 1:3
    [Ust, Ust_U, Ust_x, ~] = stagnation_state(U, x);
    v{ie} = Ust; v_u{ie} = Ust_U*reshape(dU,8,1) + Ust_x*dx;
    U = U + ep*dU; x = x + ep*dx;
  end
  check_ping(ep, v, v_u, 'stagnation state');
  
  % entire system ping
  M.oper.Re = 1e6; M.oper.alpha = 1; M.oper.Ma=0.0; rng(0);
  M.param.niglob = 10;
  solve_viscous(M); % start with a valid solution
  Nsys = size(M.glob.U,2); dU = rand(4,Nsys); dx = 0.5*rand(Nsys,1); ep = 2e-6;
  for ix = 1:2 % ping with explicit and implicit (baked-in) R_x effects
    if (ix==2), dx = dx*0; stagpoint_move(M); end % baked-in check
    v = cell(1,3); v_u = v;
    for ie = 1:3
      build_glob_sys(M); if (ix==2), jacobian_add_Rx(M); end
      v{ie} = M.glob.R; v_u{ie} = M.glob.R_U*reshape(dU,4*Nsys,1) + M.glob.R_x*dx;
      M.glob.U = M.glob.U + ep*dU; M.isol.xi = M.isol.xi + ep*dx'; 
      if (ix==2), stagpoint_move(M); end % baked-in check: stagnation point moves
    end
    M.glob.U = M.glob.U - 2*ep*dU; M.isol.xi = M.isol.xi - 2*ep*dx';       
    check_ping(ep, v, v_u, sprintf('global system, ix=%d', ix));
  end
  
  % wake system ping  
  rng(0); dU = rand(size(M.glob.U)); ep = 1e-5; v = cell(1,3); v_u = v;
  for ie = 1:3
    [R, R_U, J] = wake_sys(M, param); 
    v{ie} = R; v_u{ie} = R_U*reshape(dU(:,J),4*length(J),1);
    M.glob.U = M.glob.U + ep*dU;
  end
  M.glob.U = M.glob.U - 2*ep*dU;
  check_ping(ep, v, v_u, 'wake system');  
  
  % force calculation ping
  Nsys = size(M.glob.U,2); N = M.foil.N; rng(0); v = cell(1,3); v_u = v;
  due = rand(N,1); dU = zeros(4,Nsys); dU(4,1:N) = due'; da = 10; ep = 1e-2;
  for ie = 1:3
    calc_force(M); 
    v{ie} = M.post.cl; v_u{ie} = M.post.cl_ue*due + M.post.cl_alpha*da;
    M.glob.U = M.glob.U + ep*dU; M.oper.alpha = M.oper.alpha + ep*da;
  end
  M.glob.U = M.glob.U - 2*ep*dU; M.oper.alpha = M.oper.alpha - 2*ep*da;
  check_ping(ep, v, v_u, 'lift calculation');
  
end


%-------------------------------------------------------------------------------
% Revision history
%
% 2021-08-27: (first version) coupled solver, compressibility, cl constraint
% 2021-09-10: fixed cdf calculation for compressible flow
% 2025-04-28: added forced transition calculation
