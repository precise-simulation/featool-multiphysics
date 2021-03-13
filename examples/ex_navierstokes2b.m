function [ data, fea ] = ex_navierstokes2b( varargin )
%EX_NAVIERSTOKES2B 2D Driven cavity CFD benchmark
%
%   [ DATA, FEA ] = EX_NAVIERSTOKES2B 2D Driven cavity CFD benchmark
%
%   CFD benchmarking script comparing the OpenFOAM and SU2 CFD solvers
%   to the fully coupled FEATool Multiphysics solver, for stationary
%   laminar driven cavity test case.
%
%   See also EX_NAVIERSTOKES2, RUN_FEATOOL_BENCHMARKS

% Copyright 2013-2021 Precise Simulation, Ltd.


% Define available solvers.
solvers  = { 'FEATool', 'FEniCS', 'OpenFOAM', 'SU2' };
FEATOOL  = 1;
FENICS   = 2;
OPENFOAM = 3;
SU2      = 4;

% Define available grid cases.
grids    = { 'Quad', 'Tri', 'TriU' };
QUAD     = 1;
TRI      = 2;   % Structured triangles.
TRIU     = 3;   % Unstructured triangles.

% Define available FEM shape/basis functions.
D0       = 'sf_disc0';
D1       = 'sf_disc1';
P1       = 'sf_simp_P1';
P2       = 'sf_simp_P2';
Q1       = 'sf_quad_Q1';
Q2       = 'sf_quad_Q2';

% Set up benchmark test cases.
NLEV  = 1:5;   % Select grid levels.
cases = { FEATOOL,   TRI,   {P1,P1},  NLEV,  {};
          FEATOOL,   TRI,   {P2,P1},  NLEV(1:end-1),  {};
          FEATOOL,   QUAD,  {Q1,Q1},  NLEV,  {};
          FEATOOL,   QUAD,  {Q2,D1},  NLEV(1:end-1),  {};
          FENICS,    TRI,   {P1,P1},  NLEV(1:end-1),  {'nlrlx',0.6};
...          FENICS,    TRI,   {P2,P1},  NLEV(1:end-1),  {};
          OPENFOAM,  TRI,   {D0,D0},  NLEV,  {};
          OPENFOAM,  QUAD,  {D0,D0},  NLEV,  {};
          SU2,       TRI,   {P1,P1},  NLEV,  {};
          SU2,       QUAD,  {Q1,Q1},  NLEV,  {} };


opt.basename = mfilename();
opt.solvers  = solvers;
opt.grids    = grids;
opt.cases    = cases;
opt.fcn_fea  = @l_setup_fea_struct;
opt.fcn_err  = @l_compute_error;

[ data, fea ] = run_featool_benchmarks( opt );


%------------------------------------------------------------------------------%
function [ fea ] = l_setup_fea_struct( grid_type, sfun, i_lev )

sf_u = sfun{1};
sf_p = sfun{2};

re   = 1000;
rho  = 1;           % Density.
umax = 1;           % Maximum magnitude of inlet velocity.
l    = 1;
miu  = umax*l/re;   % Molecular/dynamic viscosity.

% Geometry definition.
gobj = gobj_rectangle( 0, l, 0, l );
gobj.boundaries = gobj.boundaries([4,1,2,3]);
fea.geom.objects = { gobj };
fea.sdim = { 'x', 'y' };

% Grid generation.
hmax = l/5;
if( grid_type==1 || grid_type==2 )
  fea.grid = rectgrid(round(l/hmax),round(l/hmax),[0 l;0 l]);
  ib = zeros(1,size(fea.grid.b,2));
  ib(fea.grid.b(3,:)==4) = 1;
  ib(fea.grid.b(3,:)==1) = 2;
  ib(fea.grid.b(3,:)==2) = 3;
  ib(fea.grid.b(3,:)==3) = 4;
  fea.grid.b(3,:) = ib;
  if( grid_type==2 )
    fea.grid = quad2tri( fea.grid );
  end
  for i=1:(i_lev-1)
    fea.grid = gridrefine( fea.grid, [] );
  end
else
  for i=1:(i_lev-1)
    hmax = hmax/2;
  end
  fea.grid = gridgen( fea, 'hmax', hmax, 'gridgen', 'gridgen2d', 'fid', [] );
end

% Problem definition.
fea = addphys(fea,@navierstokes);
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
[fea.phys.ns.sfun{1:2}] = deal( sf_u );
fea.phys.ns.sfun{3} = sf_p;

% Boundary conditions.
fea.phys.ns.bdr.sel(4) = 2;
fea.phys.ns.bdr.coef{2,end}{1,4} = umax;

fea = parsephys(fea);
fea = parseprob(fea);

%------------------------------------------------------------------------------%
function [ err ] = l_compute_error( fea )

vort = evalexpr( 'vx-uy', [0.53;0.564], fea );
err  = abs(-2.068 - vort)/2.068;
