function [ fea, out ] = ex_convdiff3( varargin )
%EX_CONVDIFF3 1D Time dependent convection and diffusion equation example.
%
%   [ FEA, OUT ] = EX_CONVDIFF3( VARARGIN ) 1D time dependendt convection and diffusion equation on
%   a line with an infinately oscillating exact solution. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       w           scalar {1.5*pi}        Simulation parameter
%       U           scalar {2}             Simulation parameter
%       a           scalar {1}             Convection velocity
%       nu          scalar {0.1}           Diffusion coefficient
%       hmax        scalar {1/25}          Max grid cell size
%       dt          scalar {0.02}          Time step size
%       ischeme     scalar {3}             Time stepping scheme
%       sfun        string {sflag2}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'w',        1.5*pi; ...
  'U',        2; ...
  'a',        1; ...
  'nu',       0.1; ...
  'hmax',     1/25; ...
  'dt'        0.02; ...
  'ischeme'   3; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      3e-2; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


w  = opt.w;
U  = opt.U;
a  = opt.a;
nu = opt.nu;
l1 = ['(',num2str(a),'+sqrt(',num2str(a^2),'+',num2str(4*nu*w),'*i))/',num2str(2*nu)];
l2 = ['(',num2str(a),'-sqrt(',num2str(a^2),'+',num2str(4*nu*w),'*i))/',num2str(2*nu)];
refsol  = ['(exp(',l1,'*x)-exp(',l2,'*x))/(exp(',l1,')-exp(',l2,'))*',num2str(U),'*exp(i*',num2str(w),'*t)'];
refsol0 = ['(exp(',l1,'*x)-exp(',l2,'*x))/(exp(',l1,')-exp(',l2,'))*',num2str(U)];


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
fea.bdr.d{1} = 0;
fea.bdr.d{2} = [num2str(U),'*cos(',num2str(w),'*t)'];
fea.bdr.n    = cell(1,2);
x = fea.grid.p';
if( strcmp( opt.sfun,'sflag2' ) )
  x = [ x; (x(2:end)+x(1:end-1))/2 ];
end
init = real( eval( refsol0 ) );
[fea.sol.u,tlist] = solvetime( fea, 'fid', fid, 'init', init, 'ischeme', opt.ischeme, 'tstep', opt.dt, 'tmax', 1 );


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
    u_r = real( eval( refsol ) );
    plot( sort(x), u_r(ix), 'r--' );
    title( ['Solution at time ',num2str(t)])
    xlabel( 'x' )
    drawnow
  end
end


% Error checking.
for i_sol=1:numel(tlist)
  u_i = fea.sol.u(:,i_sol);
  t   = tlist(i_sol);
  u_r = real( eval( refsol ) );
  errnm(i_sol) = norm( u_i - u_r )/norm( u_r );
end
out.err  = errnm;
out.pass = all( errnm<opt.tol );


if ( nargout==0 )
  clear fea out
end
