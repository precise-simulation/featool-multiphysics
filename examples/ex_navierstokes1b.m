function [ data, fea ] = ex_navierstokes1b( varargin )
%EX_NAVIERSTOKES1B 2D Channel flow CFD benchmark
%
%   [ DATA, FEA ] = EX_NAVIERSTOKES1B 2D Channel flow CFD benchmark
%
%   CFD benchmarking script comparing the OpenFOAM and SU2 CFD solvers
%   to the fully coupled and monolithic FEniCS and FEATool Multiphysics
%   solvers, for stationary laminar Poiseuille flow in a rectangular
%   channel. The inflow profile is constant and the outflow should
%   assume a parabolic profile u(y) = U_max*4/h^2*y*(h-y).
%
%   See also EX_NAVIERSTOKES1, RUN_FEATOOL_BENCHMARKS

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
cases = { FEATOOL,   TRI,   {P1,P1},  NLEV;
          FEATOOL,   TRI,   {P2,P1},  NLEV(1:end-1);
          FEATOOL,   QUAD,  {Q1,Q1},  NLEV;
          FEATOOL,   QUAD,  {Q2,D1},  NLEV(1:end-1);
          FENICS,    TRI,   {P1,P1},  NLEV;
          FENICS,    TRI,   {P2,P1},  NLEV(1:end-1);
          OPENFOAM,  TRI,   {D0,D0},  NLEV;
          OPENFOAM,  QUAD,  {D0,D0},  NLEV;
          SU2,       TRI,   {P1,P1},  NLEV;
          SU2,       QUAD,  {Q1,Q1},  NLEV };


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

rho  = 1;           % Density.
miu  = 0.001;       % Molecular/dynamic viscosity.
u_in = 2/3*0.3;     % Inlet velocity.
h    = 0.5;         % Height of rectangular domain.
l    = 25*h;        % Length of rectangular domain.

fea.refsol = ['4*',num2str(3/2*u_in),'*(y*(',num2str(h),'-y))/',num2str(h),'^2'];

% Geometry definition.
fea.geom.objects = { gobj_rectangle(0,l,0,h) };
fea.sdim = { 'x', 'y' };

% Grid generation.
hmax = min(h,l)/5;
if( grid_type==1 || grid_type==2 )
  fea.grid = rectgrid(round(l/hmax/5),round(h/hmax),[0 l;0 h]);
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
% fea.phys.ns.eqn.seqn{1} = '-miu_ns*(ux_x + uy_y) + rho_ns*(u*ux_t + v*uy_t) + p_x = 0';   % Gradient formulation.
% fea.phys.ns.eqn.seqn{2} = '-miu_ns*(vx_x + vy_y) + rho_ns*(u*vx_t + v*vy_t) + p_y = 0';
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
[fea.phys.ns.sfun{1:2}] = deal( sf_u );
fea.phys.ns.sfun{3} = sf_p;

% Boundary conditions.
fea.phys.ns.bdr.sel(4) = 2;
fea.phys.ns.bdr.sel(2) = 4;
fea.phys.ns.bdr.coef{2,end}{1,4} = u_in;

fea = parsephys(fea);
fea = parseprob(fea);

% Enforce straight-out outflow boundary (v=0).
fea.bdr.d{1}{2} = [];
fea.bdr.d{2}{2} = 0;
fea.bdr.n{1}{2} = 0;
fea.bdr.n{2}{2} = [];

%------------------------------------------------------------------------------%
function [ err ] = l_compute_error( fea )

s_err = ['abs( sqrt((',fea.refsol,')^2) - sqrt(u^2+v^2) )'];
ix = fea.grid.b(3,:) == 2;
ind_c = unique( fea.grid.b(1,ix) );   % Indices to grid cells on outflow boundary.
if( strcmp(fea.sfun{1},'sf_disc0') )
  err = evalexprc( s_err, ind_c, fea );
  ref = evalexprc( ['sqrt((',fea.refsol,')^2)'], ind_c, fea );
else
  ind_e = fea.grid.b(2,ix);
  ind_p = [ fea.grid.c( sub2ind( size(fea.grid.c), ind_e, ind_c ) ), ...
            fea.grid.c( sub2ind( size(fea.grid.c), mod(ind_e,size(fea.grid.c,1))+1, ind_c ) ) ];
  ind_p = unique(ind_p);
  err = evalexprp( s_err, fea, 1, ind_p );
  ref = evalexprp( ['sqrt((',fea.refsol,')^2)'], fea, 1, ind_p );
end
err = sqrt( sum(err.^2) / sum(ref.^2) );
