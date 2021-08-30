function [ fea, out ] = ex_convdiff6( varargin )
%EX_CONVDIFF6 1D Stationary convection and diffusion equation example.
%
%   [ FEA, OUT ] = EX_CONVDIFF6( VARARGIN ) 1D stationary convection and diffusion equation on a
%   line L with the exact solution u(x) = c0 + f/v*x + f*D/v^2*( exp(-v*L/D) - exp((x-L)*v/D) ).
%   A Dirichlet condition c0 is prescribed at x=0 and at x=L is left free dc(x=L)/dx = 0.
%   Accepts the following property/value pairs.
%
%   Reference:
%
%   [1] Case B3 in M. Th. van Genuchten and W. J. Alves, Analytical Solutions to the
%       One-Dimensional Convective-Dispersive Solution Transport Equation, pp. 29,
%       Technical Bulletin 1661, Publisher: U.S. Department of Agriculture, 1982.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       D           scalar {2}             Diffusion coefficient
%       v           scalar {3}             Convection velocity
%       f           scalar {4}             Reaction source term
%       c0          scalar {5}             Dirichlet boundary condition at x=0
%       L           scalar {6}             Length of domain
%       nx          scalar {25}            Number of grid cells
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'D',     2; ...
            'v',     3; ...
            'f',     4; ...
            'c0',    5; ...
            'L',     6; ...
            'nx',    25; ...
            'sfun',  'sflag1'; ...
            'iplot', 1; ...
            'tol',   1e-2; ...
            'fid',   1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Grid generation.
fea.grid = linegrid( opt.nx, 0, opt.L );


% Problem definition.
fea.sdim  = { 'x' };

fea = addphys( fea, @convectiondiffusion );
fea.phys.cd.sfun = { opt.sfun };

fea.phys.cd.eqn.coef{2,4} = { opt.D };
fea.phys.cd.eqn.coef{3,4} = { opt.v };
fea.phys.cd.eqn.coef{4,4} = { opt.f };

fea.phys.cd.bdr.sel(1)         = 1;
fea.phys.cd.bdr.coef{1,end}{1} = opt.c0;


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')), 'fid', opt.fid );



% Postprocessing.
refsol = sprintf( '%d + %d*x + %d*( exp(-%d) - exp((x-%d)*%d) )', ...
                  opt.c0, opt.f/opt.v, opt.f*opt.D/opt.v^2, ...
                  opt.v*opt.L/opt.D, opt.L, opt.v/opt.D );

x     = linspace( 0, opt.L,  100 );
c     = evalexpr( 'c',    x, fea );
c_ref = evalexpr( refsol, x, fea );
if( opt.iplot>0 )
  plot( x, c,     'b-.', 'linewidth', 2 );
  hold on, grid on
  plot( x, c_ref, 'r--', 'linewidth', 1 );
  legend( 'computed solution', 'reference solution', 'location', 'southeast' )
end


% Error checking.
ix = find( c_ref~=0 );
errnm    = norm( (c(ix) - c_ref(ix)) ./ c_ref(ix) );
out.err  = errnm;
out.pass = all( errnm<opt.tol );

if( nargout==0 )
  clear fea out
end
