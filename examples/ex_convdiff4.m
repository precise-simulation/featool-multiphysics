function [ fea, out ] = ex_convdiff4( varargin )
%EX_CONVDIFF4 1D Burgers equation (convection and diffusion) example.
%
%   [ FEA, OUT ] = EX_CONVDIFF4( VARARGIN ) 1D Burgers equation with steady solution,
%   u_t + (b*u-c)*u_x - nu*u_xx = 0 with exact solution c/b*(1-tanh(c/(2*nu)*(x-x0))).
%   Tests both time dependent and steady solvers. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       b           scalar {1.0}           Strength of nonlinearity
%       c           scalar {0.13}          Convection velocity
%       nu          scalar {0.01}          Diffusion coefficient
%       x0          scalar {0.5}           Posision of smooth shock
%       hmax        scalar {1/25}          Max grid cell size
%       ischeme     scalar {-1}            Solver scheme (<0 = stationary)
%       nsolve      scalar {2}             Nonlinear solver (when ischeme<0)
%       dt          scalar {0.1}           Time step size
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
            'b',        1.0; ...
            'c',        0.13; ...
            'nu',       0.01; ...
            'x0',       0.5; ...
            'hmax',     1/20; ...
            'ischeme'   -1; ...
            'nsolve'    2; ...
            'dt'        0.1; ...
            'sfun',     'sflag1'; ...
            'iplot',    1; ...
            'tol',      1e-2; ...
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


b  = opt.b;
c  = opt.c;
nu = opt.nu;
x0 = opt.x0;
refsol = [num2str(c/b),'*(1-tanh(',num2str(c/(2*nu)),'*(x-',num2str(x0),')))'];
bcd1   = c/b*(1-tanh((c/(2*nu)*(0-x0))));
bcd2   = c/b*(1-tanh((c/(2*nu)*(1-x0))));

% Grid generation.
fea.grid = linegrid( 1/opt.hmax, 0, 1 );


% Problem definition.
fea.sdim  = { 'x' };
fea = addphys( fea, @convectiondiffusion );
fea.phys.cd.sfun = { opt.sfun };
fea.phys.cd.eqn.coef{2,4} = { opt.nu };
fea.phys.cd.eqn.coef{3,4} = { [num2str(b),'*c-',num2str(c)] };
fea = parsephys(fea);


% Parse and solve problem.
fea = parseprob( fea );
fea.bdr.d{1} = bcd1;
fea.bdr.d{2} = bcd2;
fea.bdr.n    = cell(1,2);
x = fea.grid.p';
if( strcmp( opt.sfun,'sflag2' ) )
  x = [ x; (x(2:end)+x(1:end-1))/2 ];
end
if( opt.ischeme<0 )
  init = 0;
  jac.form  = {[1;1]};
  jac.coef  = {[num2str(b),'*cx']};
  fea.sol.u = solvestat( fea, 'fid', fid, 'init', init, 'maxnit', 1000, 'nsolve', 2, 'jac', jac, 'nsolve', opt.nsolve );
else
  init = [num2str(bcd1),'+x*',num2str(bcd2-bcd1)];
  [fea.sol.u,tlist] = solvetime( fea, 'fid', fid, 'init', init, 'ischeme', opt.ischeme, 'tstep', opt.dt, 'tmax', 10 );
end


% Postprocessing.
if( opt.iplot>0 )
  figure
  i_sol = size(fea.sol.u,2);
  [~,ix] = sort( x );
  postplot( fea, 'surfexpr', 'c', 'solnum', i_sol );
  hold on
  u_r = real( eval( refsol ) );
  plot( sort(x), u_r(ix), 'r--' );
  title( 'Steady solution' )
  xlabel( 'x' )
  drawnow
end


% Error checking.
i_sol = size(fea.sol.u,2);
u_i   = fea.sol.u(:,i_sol);
u_r   = real( eval( refsol ) );
errnm = norm( u_i - u_r )/norm( u_r );
out.err  = errnm;
out.pass = all( errnm<opt.tol );


if ( nargout==0 )
  clear fea out
end
