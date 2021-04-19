function [ fea, out ] = ex_heattransfer1( varargin )
%EX_HEATTRANSFER1 2D ceramic strip with radiation and convection.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER1( VARARGIN ) 2D heat transfer of a ceramic strip with
%   both radiation and convection on the top boundary.
%
%         _                q = h*(T_inf-T) + epsilon*sigma*(T_inf^4-T^4)
%         ^            +------------------+
%         |            |                  |
%         |            |                  |
%       0.01m  T=900C  |                  |  T=900C
%         |            |                  |
%         |            |                  |
%         v            +------------------+
%                           dt/dn = 0
%                      |<---- 0.02m ----->|
%
%   The ceramic has a thermal conductivity of 3 W/mC and the sides are fixed
%   at a temperature of 900C while the bottom boundary is insulated. The surrounding
%   temperature is 50C. The top boundary is exposed to both natural convection (with
%   a film coefficient h=50W/m^2K) and radiation (with emissivity epsilon=0.7
%   and the Stefan-Boltzmann 5.669e-8 W/m^2K^4). The solution is sought at three
%   points along the vertical symmetry line.
%
%   Reference:
%
%      [1] Holman, J. P., Heat Transfer, Fifth Edition, New York: McGraw-Hill,
%          1981, page 96, Example 3-8.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.001}         Grid cell size
%       igrid       scalar {0}/1/2         Cell type (0=quadrilaterals, 1=triangles,
%       solver      string fenics/{}       Use FEniCS or default solver
%       ischeme     scalar {0}             Time stepping scheme (0 = stationary)
%       sfun        string {sflag1}        Finite element shape function
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.001;
            'igrid',    0;
            'solver',   '';
            'ischeme',  0;
            'sfun',     'sflag1';
            'iplot',    1;
            'tol',      0.01;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
if( opt.ischeme==2 && ~got.tol )
  opt.tol = 0.05;
end


% Geometry definition.
gobj = gobj_rectangle( 0, 0.02, 0, 0.01 );
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case 0
    fea.grid = rectgrid( round(0.02/opt.hmax), round(0.01/opt.hmax), [0 0.02;0 0.01] );
  case 1
    fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );
  case 2
    fea.grid = rectgrid( round(0.02/opt.hmax), round(0.01/opt.hmax), [0 0.02;0 0.01] );
    fea.grid = quad2tri( fea.grid, 1 );
end


% Problem definition.
fea.sdim  = { 'x', 'y' };             % Space coordinate name.
fea = addphys( fea, @heattransfer );  % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };      % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{3,end} = 3;      % Thermal conductivity.

% Boundary conditions.
fea.phys.ht.bdr.sel = [3 1 4 1];
fea.phys.ht.bdr.coef{1,end}   = { [] 900+273 [] 900+273 };
fea.phys.ht.bdr.coef{4,end}{3}{2} = 50;
fea.phys.ht.bdr.coef{4,end}{3}{3} = 50+273;
fea.phys.ht.bdr.coef{4,end}{3}{4} = 0.7*5.669e-8;
fea.phys.ht.bdr.coef{4,end}{3}{5} = 50+273;


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', 0.1, 'tmax', 1, 'ischeme', opt.ischeme );
else
  if( opt.ischeme<=0 )
    fea.sol.u = solvestat( fea, 'fid', opt.fid, 'init', {'T0_ht'} );
  else
    [fea.sol.u,fea.sol.t] = solvetime( fea, 'fid', opt.fid, 'init', {'T0_ht'}, ...
                                       'tstep', 0.1, 'tmax', 1, 'ischeme', opt.ischeme );
  end
end

% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'T', 'isoexpr', 'T' )
  title('Temperature, T')
end


% Error checking.
T2_sol = evalexpr( 'T', [0.01;0.01], fea );
T2_ref = 984;
T5_sol = evalexpr( 'T', [0.01;0.005], fea );
T5_ref = 1064;
T8_sol = evalexpr( 'T', [0.01;0], fea );
T8_ref = 1088;
out.err  = abs([T2_sol-T2_ref T5_sol-T5_ref T8_sol-T8_ref])./[T2_ref T5_ref T8_ref];
out.pass = all(out.err<opt.tol);

if( nargout==0 )
  clear fea out
end
