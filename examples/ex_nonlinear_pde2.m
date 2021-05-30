function [ fea, out ] = ex_nonlinear_pde2( varargin )
%EX_NONLINEAR_PDE2 2D nonlinear PDE example.
%
%   [ FEA, OUT ] = EX_NONLINEAR_PDE2( VARARGIN ) Solves (x*u*u_x)_x + uy_y = 0 with
%   Dirichlet boundary conditions u(x=x1)=u1 and u(x=x2)=u2. Compares against the
%   analytical solution u(x) = sqrt(2*c1*log(x)+2*c2), with c1 = -(u1^2-u2^2)/
%   (log(x2)-log(x1))/2 and c2 = (log(x2)*u1^2-u2^2*log(x1))/(log(x2)-log(x1))/2.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {1/10}          Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       u1          scalar                 Dirichlet boundary value at x = x1
%       u2          scalar                 Dirichlet boundary value at x = x2
%       x1          scalar                 Left most domain x-coordinate
%       x2          scalar                 Right most domain x-coordinate
%       ischeme     scalar {-1}            Solver scheme (<0 = stationary)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'hmax',     1/10;
  'sfun',     'sflag1';
  'u1',       1;
  'u2',       2;
  'x1',       3;
  'x2',       5;
  'ischeme', -1;
  'iplot',    1;
  'tol',      1e-3;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Grid generation.
nx       = round( 1/opt.hmax );
fea.grid = rectgrid( nx, nx, [opt.x1 opt.x2;0 1] );

% Coefficients and expressions.
fea.sdim  = { 'x' 'y' };
fea.expr  = { 'u1'     opt.u1;
              'u2'     opt.u2;
              'x1'     opt.x1;
              'x2'     opt.x2;
              'c1'     '-(u1^2-u2^2)/(log(x2)-log(x1))/2';
              'c2'     '(log(x2)*u1^2-u2^2*log(x1))/(log(x2)-log(x1))/2';
              'u_ref'  'sqrt(2*c1*log(x)+2*c2)' };

% Physics mode and equation settings.
fea = addphys( fea, @customeqn, { 'u' } );
fea.phys.ce.sfun = { opt.sfun };
fea.phys.ce.eqn.seqn = 'u'' - x*u*ux_x - uy_y = 0';
fea.phys.ce.eqn.coef = { 'u0_ce', 'u_0', 'Initial condition for u', { 'u1+(u2-u1)*(x-x1)/(x2-x1)' } };

% Boundary settings.
fea.phys.ce.bdr.coef = { 'bcnd_ce', 'Dirichlet/Neumann boundary conditions',...
                         'Dirichlet/Neumann boundary conditions', ...
                         {'Dirichlet, r_u', 'Neumann, g_u' }, ...
                         { 0 1 0 1 }, [], { 0 'u2' 0 'u1' } };


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
if( opt.ischeme<1 )
  fea.sol.u = solvestat( fea, 'fid', fid, 'init', {'u0_ce'} );
else
  fea.sol.u = solvetime( fea, 'fid', fid, 'init', {'u0_ce'}, 'ischeme', opt.ischeme );
end


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', 'u', 'surfhexpr', 'u' )
  title('Solution u')
  u_ref = evalexprp( 'u_ref', fea );
  hold on
  plot3( fea.grid.p(1,:), fea.grid.p(2,:), u_ref, 'kx' )
end


% Error checking.
err = evalexprp( 'abs(u_ref-u)', fea );
ref = evalexprp( 'u',            fea );
err = sqrt(sum(err.^2)/sum(ref.^2));

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end

out.err  = err;
out.pass = out.err<opt.tol;
if( nargout==0 )
  clear fea out
end
