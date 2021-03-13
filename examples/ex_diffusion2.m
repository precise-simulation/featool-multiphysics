function [ fea, out ] = ex_diffusion2( varargin )
%EX_DIFFUSION2 1D Time dependent diffusion equation examples.
%
%   [ FEA, OUT ] = EX_DIFFUSION2( VARARGIN ) 1D time dependent diffusion equation on
%   a line with exact solutions. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       nu          scalar {1e-1}          Diffusion coefficient
%       icase       scalar {2}             Test case
%       hmax        scalar {1/25}          Max grid cell size
%       dt          scalar {0.1}           Time step size
%       ischeme     scalar {3}             Time stepping scheme
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'nu',       1e-1; ...
  'icase',    2; ...
  'hmax',     1/25; ...
  'dt'        0.1; ...
  'ischeme'   3; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      1e-2; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


tmax = 1;
switch( opt.icase )
  case 1
    refsol = '(1-x)*2.2+x*3.3';
  case 2
    n = 1;
    refsol = ['exp(-',num2str(opt.nu*n^2*pi^2),'*t)*sin(',num2str(n*pi),'*x)'];
end


% Grid generation.
fea.grid = linegrid( 1/opt.hmax, 0, 1 );


% Problem definition.
fea.sdim  = { 'x' };
fea = addphys( fea, @convectiondiffusion );
fea.phys.cd.sfun = { opt.sfun };
fea.phys.cd.eqn.coef{2,4} = { opt.nu };
fea = parsephys(fea);


% Parse and solve problem.
x = fea.grid.p';
n = length(x);
if( strcmp( opt.sfun,'sflag2' ) )
  x = [ x; (x(2:end)+x(1:end-1))/2 ];
end
t   = 0;
u0  = eval( refsol );
fea = parseprob( fea );
if( opt.ischeme>0 )

  fea.bdr.d{1} = refsol;
  fea.bdr.d{2} = refsol;
  fea.bdr.n    = cell(1,2);
  [fea.sol.u,tlist] = solvetime( fea, 'fid', fid, 'init', u0, 'ischeme', opt.ischeme, 'tstep', opt.dt, 'tmax', tmax, 'nstbwe', 0, 'tstop', 0, 'imass', 2 );

else

  [M,A,f] = assembleprob( fea, 'f_m', 1, 'imass', 1, 'f_a', 1, 'f_f', 1, 'f_sparse', 1 );
  M = spdiags( full(sum(M')'), 0, size(M,1), size(M,1) );
  fea.sol.u = u0;
  dt = opt.dt;
  it = 0;
  tlist = 0;
  if( opt.ischeme==-1 )       % Crank-Nicolson.
    C = [ M + dt/2*A ];
  elseif( opt.ischeme==-2 )   % Backward Euler.
    C = [ M + dt*A ];
  end
  C(1,:) = 0;
  C(1,1) = 1;
  C(n,:) = 0;
  C(n,n) = 1;
  while 1
    t  = t  + dt;
    it = it + 1;

    u_r = eval( refsol );
    if( opt.ischeme==-1 )       % Crank-Nicolson.
      b = [ [ M - dt/2*A ]*u0 + dt*f ];
    elseif( opt.ischeme==-2 )   % Backward Euler.
      b = [ M*u0 + dt*f ];
    end
    b(1) = u_r(1);
    b(n) = u_r(n);

    u1 = C\b;

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
for i_sol=1:numel(tlist)
  u_i = fea.sol.u(:,i_sol);
  t   = tlist(i_sol);
  u_r = eval( refsol );
  errnm(i_sol) = norm( u_i - u_r )/norm( u_r );
end
out.err  = errnm;
out.pass = all( errnm<opt.tol );


if ( nargout==0 )
  clear fea out
end
