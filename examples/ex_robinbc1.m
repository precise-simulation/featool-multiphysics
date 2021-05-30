function [ fea, out ] = ex_robinbc1( varargin )
%EX_ROBINBC1 1D Robin boundary condition example.
%
%   [ FEA, OUT ] = EX_ROBINBC1( VARARGIN ) Convection, diffusion,
%   and reaction equation uxx+u*ux-u=exp(2x) on a line with a Robin
%   boundary condition u(0)+ux(0)=2 and u(1)=exp(1), and exact
%   solution exp(x). Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {1/10}          Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',   1/10;
            'sfun',   'sflag1';
            'iplot',  1;
            'icub',   2;
            'refsol', 'exp(x)';
            'tol',    3e-3;
            'fid',    1 };
[got,opt] = parseopt(cOptDef,varargin{:});

nx = round( 1/opt.hmax );

fea.sdim  = {'x'};
fea.grid  = linegrid(nx,0,1);
fea.dvar  = {'u'};
fea.sfun  = {opt.sfun};
fea.eqn   = parseeqn( 'ux_x + u*ux_t - u_t = exp(2*x)', ...
                      fea.dvar, fea.sdim );
fea.bdr.d = {[],opt.refsol};
fea.bdr.n = {'2-u',[]};


% Parse and solve problem.
fea       = parseprob(fea);           % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',opt.fid,'icub',opt.icub); % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  h1 = postplot( fea, 'surfexpr', 'u', 'color', 'b' );
  h2 = postplot( fea, 'surfexpr', opt.refsol, ...
                 'linestyle', '--', 'color', 'r' );
  grid on
  legend( [h1(1),h2(1)], 'computed solution', ...
          'analytic solution', 'location', 'northwest' )
  xlabel('x')
  ylabel('u')
  axis normal
end


% Error checking.
xi = [1/2; 1/2];
s_err = ['abs(',opt.refsol,'-u)'];
err = evalexpr0(s_err,xi,1,1:size(fea.grid.c,2),[],fea);
ref = evalexpr0('u',xi,1,1:size(fea.grid.c,2),[],fea);
err = sqrt(sum(err.^2)/sum(ref.^2));

if( ~isempty(opt.fid) )
  fprintf(opt.fid,'\nL2 Error: %e\n',err)
  fprintf(opt.fid,'\n\n')
end

out.err  = err;
out.tol  = opt.tol;
out.pass = out.err<out.tol;
if( nargout==0 )
  clear fea out
end
