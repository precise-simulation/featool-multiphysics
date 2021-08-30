function [ data, fea ] = ex_navierstokes8b( varargin )
%EX_NAVIERSTOKES8B Axisymmtric laminar pipe flow benchmark
%
%   [ DATA, FEA ] = EX_NAVIERSTOKES8B Axisymmtric laminar pipe flow benchmark
%
%   Sets up and solves stationary axisymmetric laminar Poiseuille flow
%   in a circular pipe. The inflow profile is constant and the outflow
%   should assume a parabolic profile. Compares the FEATool, FEniCS,
%   and OpenFOAM solvers.
%
%   See also EX_NAVIERSTOKES8, RUN_FEATOOL_BENCHMARKS

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
          OPENFOAM,  QUAD,  {D0,D0},  NLEV };


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

rho  = 2;           % Density.
miu  = 3;           % Molecular/dynamic viscosity.
w_in = 1;           % Inlet velocity.
r    = 1;           % Radius of rectangular domain.
l    = 25*r;        % Length of rectangular domain.

fea.refsol = ['2*',num2str(w_in),'*(1-(r/',num2str(r),').^2)'];

% Geometry definition.
fea.geom.objects = { gobj_rectangle(0,r,0,l) };
fea.sdim = { 'r', 'z' };

% Grid generation.
hmax = min(r,l)/5;
if( grid_type==1 || grid_type==2 )
  fea.grid = rectgrid(round(r/hmax),round(l/hmax/5),[0 r;0 l]);
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
fea = addphys(fea,{@navierstokes,true});
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
[fea.phys.ns.sfun{1:2}] = deal( sf_u );
fea.phys.ns.sfun{3} = sf_p;

% Boundary conditions.
dtol = 1e-3;
i_in  = findbdr( fea, ['z<=',num2str(dtol)] );
i_out = findbdr( fea, ['z>=',num2str(3-dtol)] );
i_sym = findbdr( fea, ['r<=',num2str(dtol)] );
fea.phys.ns.bdr.sel(i_in) = 2;
fea.phys.ns.bdr.sel(i_out) = 4;
fea.phys.ns.bdr.coef{2,end}{2,i_in} = w_in;
fea.phys.ns.bdr.sel(i_sym) = 5;

fea = parsephys(fea);
fea = parseprob(fea);

% Enforce straight-outflow (u=0).
fea.bdr.d{1}{3} = 0;
fea.bdr.d{2}{3} = [];
fea.bdr.n{1}{3} = [];
fea.bdr.n{2}{3} = 0;

% Enforce slip on symmetry axis (u=0).
fea.bdr.d{1}{4} = 0;
fea.bdr.d{2}{4} = [];
fea.bdr.n{1}{4} = [];
fea.bdr.n{2}{4} = 0;

%------------------------------------------------------------------------------%
function [ err ] = l_compute_error( fea )

s_err = ['abs( sqrt((',fea.refsol,')^2) - sqrt(u^2+w^2) )'];
ix = fea.grid.b(3,:) == 3;
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
