function [ fea, out ] = ex_navierstokes3( varargin )
%EX_NAVIERSTOKES3 Stationary laminar 2D flow around a cylinder (Re=20).
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES3( VARARGIN ) 2D validation and CFD
%   benchmark for stationary laminar incompressible flow around a
%   cylinder (Reynolds number, Re = 20).
%
%   The drag and lift coefficients, and pressure difference between
%   the front and back of the cylinder are computed and compared with
%   reference values [1, 2].
%
%   References:
%
%   [1] John V, Matthies G. Higher-order finite element discretizations in a
%       benchmark problem for incompressible flows. International Journal for
%       Numerical Methods in Fluids, 37(8):885–903, 2001.
%
%   [2] Nabh G. On higher order methods for the stationary incompressible
%       Navier-Stokes equations. PhD Thesis, Universitaet Heidelberg, 1998.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar {2}             Grid type: >0 regular (igrid refinements)
%                                                     <0 unstruc. grid (with hmax=|igrid|)
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       solver      string {default}       Solver selection default, openfoam, su2, or fenics
%       iplot       logical false/{true}   Plot and visualize solution
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES6, EX_NAVIERSTOKES13

% Copyright 2013-2026 Precise Simulation, Ltd.

cOptDef = { 'igrid',   2;
            'sf_u',    'sflag1';
            'sf_p',    'sflag1';
            'solver',  '';
            'iplot',   true;
            'tol',     [0.05, 0.35, 0.1];
            'fid',     1;
            'iphys',   true };
[got,opt] = parseopt(cOptDef,varargin{:});


% Model parameters.
rho   = 1;      % Density.
miu   = 0.001;  % Molecular/dynamic viscosity.
umax  = 0.3;    % Maximum magnitude of inlet velocity.
umean = 0.2;    % Mean inlet velocity (2/3 umax).

cd_ref = 5.5795352338;    % Reference drag coefficient.
cl_ref = 0.010618937712;  % Reference lift coefficient.
dp_ref = 0.11752016697;   % Reference pressure difference.


% Geometry.
h    = 0.41;  % Height of rectangular domain.
l    = 2.2;   % Length of rectangular domain.
xc   = 0.2;   % x-coordinate of cylinder center.
yc   = 0.2;   % y-coordinate of cylinder center.
diam = 0.1;   % Diameter of cylinder.

fea.sdim = { 'x', 'y' };
R1 = gobj_rectangle( 0, l, 0, h, 'R1' );
C1 = gobj_circle( [xc, yc], diam/2, 'C1' );
fea.geom.objects = { R1, C1 };
fea = geom_apply_formula( fea, 'R1-C1' );


% Grid/mesh generation.
if ( opt.igrid >= 1 )
  fea.grid = l_cylbenchgrid_2d( opt.igrid );
else
  fea.grid = gridgen( fea, 'hmax', abs(opt.igrid), 'fid', opt.fid );
end
if ( strcmpi(opt.solver,'FENICS') && size(fea.grid.c,1) == 4 )
  warning( 'Converting quadrilateral mesh to triangular to support FEniCS.' )
  fea.grid = quad2tri( fea.grid );
end


% Boundary identification.
dtol = sqrt(eps) * 1e3;
ind_inflow = findbdr(fea, sprintf('x <= %.16g', dtol));  % Inflow boundary number.
s_inflow = sprintf('4*%.16g*y*(%.16g-y)/%.16g^2', umax, h, h);  % Inflow profile definition.
ind_outflow = findbdr(fea, sprintf('x >= (%.16g - %.16g)', l, dtol ));  % Outflow boundary number.
ind_cylinder = findbdr(fea, sprintf('sqrt((x-%.16g).^2+(y-%.16g).^2) <= %.16g', xc, yc, diam/2 + dtol));  % Cylinder boundary numbers.


% Problem definition.
srho = sprintf('%.16g', rho);
if ( opt.iphys )  % Use pre-defined Navier-Stokes equations physics mode.

  fea = addphys(fea, @navierstokes);
  fea.phys.ns.eqn.coef{1,end} = { rho };
  fea.phys.ns.eqn.coef{2,end} = { miu };
  fea.phys.ns.sfun = { opt.sf_u, opt.sf_u, opt.sf_p };

  % Boundary condition definitions.
  fea.phys.ns.bdr.sel(ind_inflow) = 2;
  fea.phys.ns.bdr.sel(ind_outflow) = 4;
  fea.phys.ns.bdr.coef{2,end}{1,ind_inflow} = s_inflow;
  fea = parsephys(fea);

else  % Manual equation definition.

  fea.dvar = { 'u', 'v', 'p' };  % Dependent variable names.
  fea.sfun = { opt.sf_u, opt.sf_u, opt.sf_p };  % Shape functions.

  % Define equation system.
  cvelx = [srho,'*u'];  % Convection velocity in x-direction.
  cvely = [srho,'*v'];  % Convection velocity in y-direction.
  fea.eqn.a.form = { [2 3 2 3;2 3 1 1]       [2;3]                   [1;2];
                     [3;2]                   [2 3 2 3;2 3 1 1]       [1;3];
                     [2;1]                   [3;1]                   [] };
  fea.eqn.a.coef = { {2*miu miu cvelx cvely}  miu                    -1;
                     miu                     {miu 2*miu cvelx cvely} -1;
                     1                       1                      [] };
  fea.eqn.f.form = { 1 1 1 };
  fea.eqn.f.coef = { 0 0 0 };


  % Boundary condition definitions.
  n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.
  fea.bdr.d = cell(3, n_bdr);
  [fea.bdr.d{1:2,:}] = deal(0);

  fea.bdr.d{1,ind_inflow} = s_inflow;

  [fea.bdr.d{:,ind_outflow}] = deal([]);

  fea.bdr.n = cell(3, n_bdr);
end


% Call solver.
fea = parseprob(fea);  % Check and parse problem struct.
switch ( upper(opt.solver) )

  case 'FENICS'
    fea = fenics( fea, 'fid', opt.fid, 'nproc', 1 );

  case 'OPENFOAM'
    fid = opt.fid; if( ~got.fid ), fid = []; end
    fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', opt.fid );

  case 'SU2'
    fid = opt.fid; if( ~got.fid ), fid = []; end
    fea.sol.u = su2( fea, 'fid', fid, 'logfid', opt.fid );


  otherwise  % Default
    jacobian.form = {[1;1] [1;1] []; [1;1] [1;1] []; [] [] []};  % Analytic Newton jacobian.
    jacobian.coef = {[srho,'*ux'] [srho,'*uy'] []; [srho,'*vx'] [srho,'*vy'] []; [] [] []};
    nlrlx = '(1 + (it > 2)) / 2';  % Dynamic non-linear relaxation parameter (it = iteration number).
    solve_param = {'nlrlx', nlrlx, 'nsolve', 2, 'jac', jacobian, 'linsolv', ''};

    fea.sol.u = solvestat( fea, 'fid', opt.fid, solve_param{:} );
end


% Visualization.
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot )
  figure
  subplot(2,1,1)
  postplot( fea, 'surfexpr', s_velm )
  title( 'Velocity field' )
  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'p' )
  title( 'Pressure' )
end


% Postprocessing.
[cd_l, cd_v, cl_l, cl_v, dp] = calc_benc_quants( fea, opt, ind_cylinder );

if ( ~isempty(opt.fid) )
  fprintf(opt.fid,'\n\nBenchmark quantities:\n\n')

  fprintf(opt.fid,'Drag coefficient,     cd = %6f (l), %6f (v) (Reference: %6f)\n', cd_l, cd_v, cd_ref)
  fprintf(opt.fid,'Lift coefficient,     cl = %6f (l), %6f (v) (Reference: %6f)\n', cl_l, cl_v, cl_ref)
  fprintf(opt.fid,'Pressure difference,  dp = %6f                   (Reference: %6f)\n', dp, dp_ref)
end

% Error checking.
out.cd = [cd_l, cd_v];
out.cl = [cl_l, cl_v];
out.dp = dp;
out.err = [abs(out.cd - cd_ref) / cd_ref;
           abs(out.cl - cl_ref) / cl_ref;
           nan, abs(out.dp - dp_ref) / dp_ref];
out.pass = all(out.err(:,2) < opt.tol(:));


if ( ~nargout )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ cd_l, cd_v, cl_l, cl_v, dp ] = calc_benc_quants( fea, opt, ind_cylinder )

rho = '1';
miu = '0.001';
umean = '0.2';
diam = '0.1';

% Calculate benchmark quantities (line integration method).
s_tfx = ['nx*p+',miu,'*(-2*nx*ux-ny*(uy+vx))'];
s_tfy = ['ny*p+',miu,'*(-nx*(vx+uy)-2*ny*vy)'];
s_cd = ['2*(',s_tfx,')/(',rho,'*',umean,'^2*',diam,')'];
s_cl = ['2*(',s_tfy,')/(',rho,'*',umean,'^2*',diam,')'];

cd_l = intbdr(s_cd, fea, ind_cylinder, 10);
cl_l = intbdr(s_cl, fea, ind_cylinder, 10);
dp = evalexpr('p', [0.15 0.25; 0.2 0.2], fea);
dp = dp(1) - dp(2);


% Calculate benchmark quantities (volume integration method).
bdrm = fea.bdr.bdrm{1};
ind_b = [];
ind_bm = [];
for ii=ind_cylinder
  ind_b  = [ind_b,  find(fea.grid.b(3,:) == ii)];
  ind_bm = [ind_bm, find(bdrm(3,:) == ii)];
end
ind_c = fea.grid.b(1, ind_b);
ind_gdof = bdrm(4, ind_bm);

% Create field 'a' with values one on the cylinder and zero everywhere else.
fea.dvar = [ fea.dvar, {'a'}       ];
fea.sfun = [ fea.sfun, fea.sfun(1) ];
fea = parseprob(fea);
n_dof = max(fea.eqn.dofm{1}(:));
u_a  = zeros(n_dof, 1);
u_a(ind_gdof) = 1;
fea.sol.u = [fea.sol.u; u_a];
fea.eqn = struct;
fea.bdr = struct;
fea = parseprob(fea);

s_tfx = ['ax*p+',miu,'*(-2*ax*ux-ay*(uy+vx))-(u*ux+v*uy)*a'];
s_tfy = ['ay*p+',miu,'*(-ax*(vx+uy)-2*ay*vy)-(u*vx+v*vy)*a'];
s_cd = ['2*(',s_tfx,')/(',rho,'*',umean,'^2*',diam,')'];
s_cl = ['2*(',s_tfy,')/(',rho,'*',umean,'^2*',diam,')'];
cd_v = intsubd(s_cd, fea, [], [], 3);
cl_v = intsubd(s_cl, fea, [], [], 3);

%------------------------------------------------------------------------------%
function [ grid ] = l_cylbenchgrid_2d( nlev )
% CYLBENCHGRID_2D Generate quadrilateral grid for the 2D DFG cylinder benchmark.

ns = 8 * 2^(nlev-1);
r  = [0.05, 0.06, 0.08, 0.11, 0.15];
x  = [0.41, 0.5, 0.7, 1, 1.4, 1.8, 2.2];
for ilev=2:nlev
  r = sort( [ r (r(1:end-1)+r(2:end))/2 ] );
  x = sort( [ x (x(1:end-1)+x(2:end))/2 ] );
end

grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
grid2 = holegrid( ns, 2^(nlev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
grid3 = rectgrid( x, ns, [0.41, 2.2; 0, 0.41] );
grid  = gridmerge( grid3, 4, grid2, 6 );

grid.s(:) = 1;
