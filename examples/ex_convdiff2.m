function [ fea, out ] = ex_convdiff2( varargin )
%EX_CONVDIFF2 1D Time dependent convection and diffusion equation example.
%
%   [ FEA, OUT ] = EX_CONVDIFF2( VARARGIN ) 1D time dependent convection and diffusion equation on
%   a line with exact solution exp(-k^2*nu*t)*sin(k*(x-a*t)) and periodic boundary conditions.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       k           scalar {2*pi}          Simulation parameter
%       a           scalar {1}             Convection velocity
%       nu          scalar {0.1}           Diffusion coefficient
%       hmax        scalar {1/25}          Max grid cell size
%       dt          scalar {0.01}          Time step size
%       ischeme     scalar {2}             Time stepping scheme
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'k',        2*pi; ...
  'a',        1; ...
  'nu',       0.1; ...
  'hmax',     1/25; ...
  'dt'        0.01; ...
  'ischeme'   2; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      1e-1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


refsol = ['exp(-',num2str(opt.k^2*opt.nu),'*t)*sin(',num2str(opt.k),'*(x-',num2str(opt.a),'*t))'];


% Grid generation.
fea.grid = linegrid( 1/opt.hmax, 0, 1 );


% Problem definition.
fea.sdim  = { 'x' };
fea = addphys( fea, @convectiondiffusion );
fea.phys.cd.sfun = { opt.sfun };
fea.phys.cd.eqn.coef{2,4} = { opt.nu };
fea.phys.cd.eqn.coef{3,4} = { opt.a  };
fea = parsephys(fea);


% Parse and solve problem.
fea = parseprob( fea );


x = fea.grid.p';
n = length(x);
if( strcmp( opt.sfun,'sflag2' ) )
  x = [ x; (x(2:end)+x(1:end-1))/2 ];
end
t  = 0;
u0 = eval( refsol );


% Assembly.
[M,A,f] = assembleprob( fea, 'f_m', 1, 'imass', 1, 'f_a', 1, 'f_f', 1, 'f_sparse', 1 );
M = spdiags( full(sum(M')'), 0, size(M,1), size(M,1) );
fea.sol.u = u0;
dt = opt.dt;
it = 0;
tlist = 0;
tmax = 1;
if( opt.ischeme==2 )       % Crank-Nicolson.
  C = [ M + dt/2*A ];
elseif( opt.ischeme==1 )   % Backward Euler.
  C = [ M + dt*A ];
end
v0 = zeros(size(C,1),1);
v0(1) =  1;
v0(n) = -1;
C = [C v0; v0' 0];


% Solver loop.
while 1
  t  = t  + dt;
  it = it + 1;

  u_r = eval( refsol );
  if( opt.ischeme==2 )       % Crank-Nicolson.
    b = [ [ M - dt/2*A ]*u0 + dt*f ];
  elseif( opt.ischeme==1 )   % Backward Euler.
    b = [ M*u0 + dt*f ];
  end

  u1 = C\[b;0];
  u1(end) = [];

  if( t>=tmax )
    break
  end

  err = norm( u1 - u_r )/norm( u_r );
  errnm(it) = err;
  if( ~isempty(fid) )
    fprintf( fid, 'Time = %f, error norm = %d\n', t, err );
  end

  fea.sol.u = [ fea.sol.u u1 ];
  tlist = [ tlist t ];
  u0 = u1;
end


% Postprocessing.
if( opt.iplot>0 )
  figure;
  if( opt.iplot>1 )
    i_sol_list = 1:numel(tlist);
  else
    i_sol_list = numel(tlist);
  end
  [~,ix] = sort( x );
  for i_sol=i_sol_list
    t = tlist(i_sol);
    clf
    postplot( fea, 'surfexpr', 'c', 'solnum', i_sol );
    hold on
    u_r = eval( refsol );
    plot( sort(x), u_r(ix), 'r--' );
    title( ['Solution at time ',num2str(t)])
    xlabel( 'x' )
    drawnow
  end
end


% Error checking.
out.err  = errnm;
out.pass = all( errnm<opt.tol );


if ( nargout==0 )
  clear fea out
end
