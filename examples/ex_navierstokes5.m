function [ fea, out ] = ex_navierstokes5( varargin )
%EX_NAVIERSTOKES5 2D Vortex flow with analytical solution.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES5( VARARGIN ) Sets up and solves time dependent
%   flow of a decaying vortex for which an analytical solution is known.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       nu          scalar {0.01}          Kinematic viscosity
%       k           scalar {pi/2}          Vortex number
%       xdim        {[-1 3]}               Domain min/max x-coordinates
%       ydim        {[-1 3]}               Domain min/max y-coordinates
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.2}           Max grid cell size
%       dt          scalar {0.01}          Time step size
%       tmax        icalar {0.1}           Simluation duration
%       ischeme     scalar {3}             Time stepping scheme
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       solver      string 'openfoam'/{''} Use OpenFOAM or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'nu',       0.01;
  'k',        pi/2;
  'xdim',     [-1 3];
  'ydim',     [-1 3];
  'igrid',    0;
  'hmax',     0.2;
  'dt',       0.01;
  'tmax',     0.1;
  'ischeme',  2;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'solver',   '';
  'iplot',    1;
  'tol',      [0.03, 0.03, 0.12];
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Exact solutions.
k = num2str(opt.k);
u_ref = ['-cos(',k,'*x).*sin(',k,'*y)*exp(-',num2str(2*opt.nu*opt.k^2),'*t)'];
v_ref = [' sin(',k,'*x).*cos(',k,'*y)*exp(-',num2str(2*opt.nu*opt.k^2),'*t)'];
p_ref = ['-1/4*(cos(2*',k,'*x)+cos(2*',k,'*y))*exp(-4*',num2str(opt.nu*opt.k^2),'*t)'];


% Grid generation.
fea.sdim = { 'x' 'y' };
if ( opt.igrid==0 )
  fea.grid = rectgrid( round(abs(diff(opt.xdim))/opt.hmax), round(diff(abs(opt.ydim))/opt.hmax), [opt.xdim;opt.ydim] );
else
  fea.geom.objects = { gobj_rectangle( opt.xdim(1), opt.xdim(2), opt.ydim(1), opt.ydim(2) ) };
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
end


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { 1 };
fea.phys.ns.eqn.coef{2,end} = { opt.nu };
fea.phys.ns.eqn.coef{5,end} = { u_ref };
fea.phys.ns.eqn.coef{6,end} = { v_ref };
fea.phys.ns.eqn.coef{7,end} = { p_ref };
fea.phys.ns.bdr.sel = [2 2 2 2];
[fea.phys.ns.bdr.coef{2,end}{1,:}] = deal( u_ref );
[fea.phys.ns.bdr.coef{2,end}{2,:}] = deal( v_ref );
[fea.phys.ns.bdr.coef{2,end}{3,:}] = deal( p_ref );
fea.phys.ns.sfun = { opt.sf_u opt.sf_u opt.sf_p };
fea = parsephys( fea );

% Remove integratl constraints.
if( isfield(fea,'constr') )
  fea = rmfield(fea,'constr');
end

% Parse and solve problem.
fea = parseprob( fea );
[fea.bdr.d{1}{:}]  = deal( u_ref );   % Exact boundary conditions.
[fea.bdr.d{2}{:}]  = deal( v_ref );
[fea.bdr.d{3}{:}]  = deal( p_ref );
fea.pnt = struct;
fea.bdr.n = cell(3,max(fea.grid.b(3,:)));
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  [fea.sol.u,fea.sol.t] = openfoam( fea, 'fid', fid, 'logfid', logfid, 'ddtScheme', 'CrankNicolson', 'deltaT', opt.dt, 'endTime', opt.tmax );
  fid = logfid;
else
  [fea.sol.u,fea.sol.t] = solvetime( fea, 'fid', fid, 'init', { u_ref v_ref p_ref }, 'ischeme', opt.ischeme, 'tstep', opt.dt, 'tmax', opt.tmax, 'nstbwe', 0, 'icub', 2 );
end


% Postprocessing.
if ( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u', 'v'}, 'arrowcolor', 'k' )
  title('Velocity field')

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'p' )
  title('Pressure')
end


% Error checking.
x = fea.grid.p(1,:)';
y = fea.grid.p(2,:)';
n = length(fea.sol.t);
for i_sol=1:n
  t = fea.sol.t(i_sol);

  u_i = evalexpr( 'u', fea.grid.p, fea, i_sol );
  v_i = evalexpr( 'v', fea.grid.p, fea, i_sol );
  p_i = evalexpr( 'p', fea.grid.p, fea, i_sol );

  u_r = eval( u_ref );
  v_r = eval( v_ref );
  p_r = eval( p_ref );

  errnm(i_sol,1) = norm( u_i - u_r )/norm( u_r );
  errnm(i_sol,2) = norm( v_i - v_r )/norm( v_r );
  errnm(i_sol,3) = norm( p_i - p_r )/norm( p_r );
end
out.err  = errnm;
out.pass = all( errnm(:) < reshape(repmat(opt.tol,n,1),[],1) );


if ( nargout==0 )
  clear fea out
end
