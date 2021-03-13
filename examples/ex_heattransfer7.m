function [ fea, out ] = ex_heattransfer7( varargin )
%EX_HEATTRANSFER7 1D Transient heat diffusion with analytic solution.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER7( VARARGIN ) Transient heat
%   diffusion problem with analytic solution. A 1 m rod is kept at
%   fixed temperature on one end and constant outward heat flux at the
%   other end as in the following illustration.
%
%              +---------- L=1m ----------+ T = 25
%            q_n = 1       T(t=0) = 25
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.1}           Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       ischeme     scalar {2}/1/3         Time stepping scheme
%       tmax        scalar {0.2}           Maximum time
%       tstep       scalar {0.01}          Time step size
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_HEATTRANSFER8.

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.1;
            'sfun',     'sflag1';
            'solver',   '';
            'ischeme',  2;
            'tmax',     0.2;
            'tstep',    0.01;
            'iplot',    1;
            'tol',      1e-3;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Grid generation.
fea.grid = linegrid( round(1/opt.hmax), 0, 1 );


% Problem definition.
fea.sdim  = { 'x' };                      % Space coordinate name.
fea = addphys( fea, @heattransfer );      % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };          % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{1,end} = 1;          % Density.
fea.phys.ht.eqn.coef{2,end} = 1;          % Heat capacity.
fea.phys.ht.eqn.coef{3,end} = 1;          % Thermal conductivity.
fea.phys.ht.eqn.coef{6,end} = { 25 };     % Initial temperature.

% Boundary conditions.
fea.phys.ht.bdr.sel = [ 4 1 ];
fea.phys.ht.bdr.coef{1,end} = { [] 25 };
fea.phys.ht.bdr.coef{4,end}{1}{1} = -1;


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', opt.tstep, 'tmax', opt.tmax, 'ischeme', opt.ischeme );
  tlist = fea.sol.t;
else
  [fea.sol.u, tlist] = solvetime( fea, 'fid', opt.fid, 'init', {'T0_ht'}, 'ischeme', opt.ischeme, ...
                                  'tmax', opt.tmax, 'tstep', opt.tstep );
end

% Postprocessing.
T_ref = refsol( fea.grid.p', tlist(end) );
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'T', 'axequal', 0 )
  title(['Temperature at t=',num2str(tlist(end))])
  xlabel('x')
  ylabel('T')

  hold on
  plot( fea.grid.p, T_ref, 'r--' )
end


% Error checking.
T_sol = evalexpr( 'T', fea.grid.p, fea );
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
  u   = u0 + 8/(1-2*n)^2/pi^2*cos((n-1/2)*pi*x)*exp(-((n-1/2)^2*pi^2)*t);
  bdo = any( max( u - u0 ) > dtol );
  u0  = u;

  if( n>1e3 )
    warning( ['Reference solution did not converge to tolerance ',num2str(dtol)] )
    break
  end
end
