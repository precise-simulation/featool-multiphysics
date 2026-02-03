function [ fea, out ] = ex_navierstokes13( varargin )
%EX_NAVIERSTOKES13 Stationary laminar 3D flow around a cylinder (Re=20).
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES13( VARARGIN ) 3D validation and CFD
%   benchmark for stationary laminar incompressible flow around a
%   cylinder (Reynolds number, Re = 20).
%
%   The drag and lift coefficients, and pressure difference between
%   the front and back of the cylinder are computed and compared with
%   reference values [1].
%
%   Reference:
%
%   [1] Braack M., Richter T. Solutions of 3D Navier-Stokes benchmark
%       problems with adaptive finite elements. Computers and Fluids,
%       35(4):372–392, 2006.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar {2}             Grid type: >0 regular (igrid refinements)
%                                                     <0 unstruc. grid (with hmax=|igrid|)
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       solver      string {default}       Solver selection default, openfoam, or su2
%       iplot       logical false/{true}   Plot and visualize solution
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES3, EX_NAVIERSTOKES6

% Copyright 2013-2026 Precise Simulation, Ltd.

cOptDef = { 'igrid',     2;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'solver',   '';
            'iplot',    true;
            'tol',      [0.2, 6, 0.25];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid = opt.fid;


% Model parameters.
rho   = 1;      % Density.
miu   = 0.001;  % Molecular/dynamic viscosity.
umax  = 0.45;   % Maximum magnitude of inlet velocity.

cd_ref = 6.185333;  % Reference drag coefficient.
cl_ref = 0.009401;  % Reference lift coefficient.
dp_ref = 0.171007;  % Reference pressure difference.


% Geometry.
h    = 0.41;  % Height of rectangular domain.
l    = 2.2;   % Length of rectangular domain.
xc   = 0.5;   % x-coordinate of cylinder center.
zc   = 0.2;   % z-coordinate of cylinder center.
diam = 0.1;   % Diameter of cylinder.

fea.sdim = { 'x', 'y', 'z' };
B1 = gobj_block( 0, l, 0, h, 0, h, 'B1' );
C1 = gobj_cylinder( [xc, 0, zc], diam/2, h, 2, 'C1' );
fea.geom.objects = { B1, C1 };
fea = geom_apply_formula( fea, 'B1-C1' );


% Grid/mesh generation.
if ( opt.igrid >= 1 )
  fea.grid = l_cylbenchgrid_3d( opt.igrid );
else
  fea.grid = gridgen( fea, 'hmax', abs(opt.igrid), 'fid', opt.fid );
end


% Problem definition.
srho = sprintf('%.16g', rho);
smiu = sprintf('%.16g', miu);
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
fea.phys.ns.sfun = { opt.sf_u, opt.sf_u, opt.sf_u, opt.sf_p };


% Boundary condition definitions.
fea.phys.ns.bdr.sel(5) = 2;
fea.phys.ns.bdr.coef{2,end}{1,5} = sprintf('16*%.16g*(y*z*(0.41-y)*(0.41-z))/0.41^4', umax);
fea.phys.ns.bdr.sel(3) = 4;
% fea.phys.ns.prop.artstab.iupw = 4;


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
switch ( upper(opt.solver) )

  case 'OPENFOAM'
    fid = opt.fid; if( ~got.fid ), fid = []; end
    fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', opt.fid );

  case 'SU2'
    fid = opt.fid; if( ~got.fid ), fid = []; end
    fea.sol.u = su2( fea, 'fid', fid, 'logfid', opt.fid );

  otherwise
    jacobian.form = {[1;1] [1;1] []; [1;1] [1;1] []; [] [] []};
    jacobian.coef = {[srho,'*ux'] [srho,'*uy'] [srho,'*uz'] [];
                     [srho,'*vx'] [srho,'*vy'] [srho,'*vz'] [];
                     [srho,'*wx'] [srho,'*wy'] [srho,'*wz'] [];
                     [] [] [] []};
    nlrlx = '(1 + (it > 2)) / 2';  % Dynamic non-linear relaxation parameter (it = iteration number).
    solve_param = {'nlrlx', nlrlx, 'nsolve', 2, 'jac', jacobian, 'linsolv', ''};

    fea.sol.u = solvestat( fea, 'fid', opt.fid, solve_param{:} );
end


% Visualization.
if ( opt.iplot )
  postplot( fea, 'sliceexpr', 'sqrt(u^2+v^2+w^2)' )
end


% Postprocessing.
s_tfx = ['nx*p+',smiu,'*(-2*nx*ux-nz*(uz+vx))'];
s_tfz = ['ny*p+',smiu,'*(-nx*(vx+uz)-2*nz*vz)'];
s_cd  = ['2*(',s_tfx,')/(',srho,'*0.2^2*0.1*0.41)'];
s_cl  = ['2*(',s_tfz,')/(',srho,'*0.2^2*0.1*0.41)'];
cd_b  = intbdr(s_cd, fea, [7:10], 10);
cl_b  = intbdr(s_cl, fea, [7:10], 10);
dp    = evalexpr('p',[0.45, 0.55; 0.205 0.205;0.21 0.21], fea);
dp = dp(1) - dp(2);

if ( ~isempty(opt.fid) )
  fprintf(opt.fid,'\n\nBenchmark quantities:\n\n')

  fprintf(opt.fid,'Drag coefficient,     cd = %6f  (Reference: %6f)\n', cd_b, cd_ref)
  fprintf(opt.fid,'Lift coefficient,     cl = %6f  (Reference: %6f)\n', cl_b, cl_ref)
  fprintf(opt.fid,'Pressure difference,  dp = %6f  (Reference: %6f)\n', dp, dp_ref)
end

% Error checking.
out.cd = cd_b;
out.cl = cl_b;
out.dp = dp;
out.err = [abs(out.cd - cd_ref) / cd_ref;
           abs(out.cl - cl_ref) / cl_ref;
           abs(out.dp - dp_ref) / dp_ref];
out.pass = all(out.err < opt.tol(:));


if ( ~nargout )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ grid ] = l_cylbenchgrid_3d( nlev )
% CYLBENCHGRID_3D Generate quadrilateral grid for the 3D DFG cylinder benchmark.

ns = 8 * 2^(nlev-1);
r  = [0.05, 0.06, 0.08, 0.11, 0.15];
x  = [0.41, 0.5, 0.7, 1, 1.4, 1.8, 2.2] + 0.3;
for ilev=2:nlev
  r = sort( [ r (r(1:end-1)+r(2:end))/2 ] );
  x = sort( [ x (x(1:end-1)+x(2:end))/2 ] );
end

grid1 = ringgrid( r, 4*ns, [], [], [0.5;0.2] );
grid2 = holegrid( ns, 2^(nlev-1), [0.3 0.71;0 0.41], 0.15, [0.5;0.2] );
grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
grid3 = rectgrid( x, ns, [0.71, 2.5; 0, 0.41] );
grid = gridmerge( grid3, 4, grid2, 6 );

grid4 = rectgrid( 1, ns, [0, 0.3; 0, 0.41] );
grid = gridmerge( grid, 10, grid4, 2 );
grid = gridextrude( grid, 5, 0.41, 2 );
grid.p(3,:) = grid.p(3,:) + 0.41;

% Renumber boundaries to match geometry.
grid.b(3,grid.b(3,:) == 6) = -7;
grid.b(3,grid.b(3,:) == 5) = -8;
grid.b(3,grid.b(3,:) == 4) = -9;
grid.b(3,grid.b(3,:) == 7) = -10;
[~,ind] = findbdr( grid, 'x<=sqrt(eps)', false );
grid.b(3,ind) = 5;
[~,ind] = findbdr( grid, 'x>=2.5-sqrt(eps)', false );
grid.b(3,ind) = 3;
[~,ind] = findbdr( grid, 'y<=sqrt(eps)', false );
grid.b(3,ind) = 2;
[~,ind] = findbdr( grid, 'y>=0.41-sqrt(eps)', false );
grid.b(3,ind) = 4;
[~,ind] = findbdr( grid, 'z<=sqrt(eps)', false );
grid.b(3,ind) = 1;
[~,ind] = findbdr( grid, 'z>=0.41-sqrt(eps)', false );
grid.b(3,ind) = 6;
grid.b(3,grid.b(3,:) == -7)  = 7;
grid.b(3,grid.b(3,:) == -8)  = 8;
grid.b(3,grid.b(3,:) == -9)  = 9;
grid.b(3,grid.b(3,:) == -10) = 10;

grid.s(:) = 1;
