function [ fea, out ] = ex_custom_equation1( varargin )
%EX_CUSTOM_EQUATION1  1D Black-Scholes custom equation example.
%
%   [ FEA, OUT ] = EX_CUSTOM_EQUATION1( VARARGIN ) 1D Black-Scholes model equation example
%   using the custom equation physics mode. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       icase       scalar {1}/2           Test case equation to solve
%       tmax        scalar {1}             Maximum/stopping time
%       len         scalar {1}             Length of domain
%       hmax        scalar {1/20}          Grid cell size
%       ischeme     scalar {3}             Time stepping scheme
%       sfun        string {sflag1}        Finite element shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'icase',    1; ...
  'tmax',     1; ...
  'len',      1; ...
  'hmax',     1/20; ...
  'ischeme'   3; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'dvname',   'u'; ...
  'tol',      1.1e-2; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;

% Grid generation.
nx       = round( opt.len/opt.hmax );
fea.grid = linegrid( nx, 0, opt.len );


% Problem definition.
u = opt.dvname;
switch( opt.icase )
  case 1
    seqn = [u,''' - 1/2*',u,'x_x - ',u,'x_t + ',u,'_t = (x-t)^5 - 10*(x-t)^4 - 10*(x-t)^3'];
  case 2
    seqn = [u,''' - 1/2*x^2*',u,'x_x - x*',u,'x_t + ',u,'_t = (x-t)^5 - 5*(x-t)^4 - 5*x*(x-t)^4 - 10*x^2*(x-t)^3'];
end
refsol = '(x-t).^5';
init_u = 'x^5';


% Set up problem struct.
fea.sdim  = { 'x' };
fea = addphys( fea, @customeqn );
fea.phys.ce.dvar = { u };
fea.phys.ce.eqn.seqn = seqn;
fea.phys.ce.sfun = { opt.sfun };


% Dirichlet BCs for the left and right boundaries.
fea.phys.ce.bdr.coef{1,5} = { 1 1 };
fea.phys.ce.bdr.coef{1,7} = { '-t^5' ['(',num2str(opt.len),'-t)^5'] };


% Check and parse problem struct.
fea = parsephys( fea );
fea = parseprob( fea );


% Call to time-dependent solver.
[fea.sol.u,tlist] = solvetime( fea, 'fid',     fid, ...
                                    'tmax',    opt.tmax, ...
                                    'init',    init_u, ...
                                    'icub',    6, ...
                                    'ischeme', opt.ischeme, ...
                                    'tstep',   opt.tmax/100);

% Postprocessing.
fea.sol.t = tlist(end);
fea.sol.u = fea.sol.u(:,end);
refsol    = strrep( refsol, 't', num2str(fea.sol.t) );
if( opt.iplot>0 )

  figure
  postplot( fea, 'surfexpr', u, 'axequal', 'off', 'linewidth', 2 )
  title( ['Solution at time ',num2str(fea.sol.t)] )
  hold on
  x = linspace(0,opt.len,25);
  u_ref = eval(refsol);
  plot( x, u_ref, 'r--' )
end


% Error checking.
err = evalexpr( ['abs(',refsol,'-',u,')'], linspace(0,opt.len,10), fea );
err = norm(err);

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end

out.err  = err;
out.pass = out.err<opt.tol;
if ( nargout==0 )
  clear fea out
end
