function [ fea, out ] = ex_heattransfer8( varargin )
%EX_HEATTRANSFER8 Space-time heat diffusion with analytic solution.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER8( VARARGIN ) One dimensional
%   transient heat diffusion problem converted to a 2D space-time
%   finite element formulation with analytic solution. A 1 m rod is
%   kept at fixed temperature on one end and constant outward heat
%   flux at the other end as in the following illustration.
%
%              +---------- L=1m ----------+ T = 25
%            q_n = 1       T(t=0) = 25
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.1}           Grid cell size in x-direction
%       igrid       scalar {0}/1/2         Cell type (0=quadrilaterals, 1=triangles,
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       tmax        scalar {0.2}           Maximum time
%       tstep       scalar {0.01}          Time step discretization size (y-direction)
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_HEATTRANSFER7.

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.025;
            'igrid',    0;
            'sfun',     'sflag1';
            'solver',   '';
            'tmax',     0.2;
            'tstep',    0.01;
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry definition.
gobj = gobj_rectangle( 0, 1, 0, opt.tmax );
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case 0
    fea.grid = rectgrid( round(1/opt.hmax), ceil(opt.tmax/opt.tstep), [0, 1; 0, opt.tmax] );
  case 1
    fea.grid = gridgen( fea, 'hmax', min(opt.hmax,opt.tstep), 'fid', opt.fid );
  case 2
    fea.grid = rectgrid( round(1/opt.hmax), ceil(opt.tmax/opt.tstep), [0, 1; 0, opt.tmax] );
    fea.grid = quad2tri( fea.grid, 1 );
end


% Problem definition.
fea.sdim  = { 'x', 'yt' };                 % Space and time coordinate names.
fea = addphys( fea, @heattransfer );      % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };          % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{1,end} = 1;          % Density (rho_ht).
fea.phys.ht.eqn.coef{2,end} = 1;          % Heat capacity (cp_ht).
fea.phys.ht.eqn.coef{3,end} = 1;          % Thermal conductivity (k_ht).
fea.phys.ht.eqn.coef{4,end} = 0;          % Convection in x-direction (u_ht).
fea.phys.ht.eqn.coef{5,end} = 1;          % Time convection coefficient (v_ht).
fea.phys.ht.eqn.coef{6,end} = 0;          % Heat source term (q_ht).
fea.phys.ht.eqn.coef{7,end} = { 25 };     % Initial temperature.

% Redefine equation.
fea.phys.ht.eqn.seqn = '-k_ht*Tx_x + rho_ht*cp_ht*u_ht*Tx_t + rho_ht*cp_ht*v_ht*Tyt_t = q_ht';

% Boundary conditions.
fea.phys.ht.bdr.sel = [ 1 1 2 4 ];
fea.phys.ht.bdr.coef{1,end} = { 25 25 [] [] };
fea.phys.ht.bdr.coef{4,end}{4}{1} = -1;

% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid );
else
  fea.sol.u = solvestat( fea, 'fid', opt.fid );
end

% Postprocessing.
T_ref = refsol( fea.grid.p(1,:)', fea.grid.p(2,:)' );
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'T', 'surfhexpr', 'T', 'boundary', 'off' )
  view(3)
  title( 'Computed Temperature' )
  xlabel('x')
  ylabel('time')
  zlabel('T')
  axis( [0 1 0 opt.tmax 0.95*min(fea.sol.u) 1.05*max(fea.sol.u)] )
  axis tight
  grid on

  subplot(1,2,2)
  fea.vars(1).name  = 'T_ref';
  fea.vars(1).descr = 'Reference Temperature';
  fea.vars(1).data  = T_ref(:);
  postplot( fea, 'surfexpr', 'T_ref', 'surfhexpr', 'T_ref', 'boundary', 'off' )
  view(3)
  title( 'Reference Temperature' )
  xlabel('x')
  ylabel('time')
  zlabel('T_{ref}')
  axis( [0 1 0 opt.tmax 0.95*min(fea.sol.u) 1.05*max(fea.sol.u)] )
  axis tight
  grid on

  rotate3d('on')
end


% Error checking.
T_sol = evalexprp( 'T', fea );
out.err  = norm( abs(T_sol-T_ref)/T_ref );
out.pass = out.err<opt.tol;

if( nargout==0 )
  clear fea out
end


% -----------------------------------
function [ u ] = refsol( x, t, dtol )

if( nargin<3 )
  dtol = 1e-7;
end

u0  = x + 24;
bdo = true;
n   = 0;
while( bdo )

  n   = n + 1;
  u   = u0 + 8/(1-2*n)^2/pi^2*cos((n-1/2)*pi*x).*exp(-((n-1/2)^2*pi^2)*t);
  bdo = any( max( u - u0 ) > dtol );
  u0  = u;

  if( n>2e3 )
    warning( ['Reference solution did not converge to tolerance ',num2str(dtol)] )
    break
  end
end
