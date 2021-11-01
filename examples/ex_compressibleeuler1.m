function [ fea, out ] = ex_compressibleeuler1( varargin )
%EX_COMPRESSIBLEEULER1 SOD Shock tube example.
%
%   [ FEA, OUT ] = EX_COMPRESSIBLEEULER1( VARARGIN ) Sets up and
%   solves an instationary SOD shock tube compressible Euler equation
%   example on the unit line. Initial data are rhoL, uL=0, pL=1 and
%   rhoR=0.125, uR=0, pR=0.1 with an inital discontinuity at t=0 and
%   x=0.5. The computed solution with shocks, expansion, and compression
%   waves is compared with the analytical Riemann solution.
%
%   Ref. G.A. Sod, A Survey of Several Finite Difference Methods for
%   Systems of Nonlinear Hyperbolic Conservation Laws, JCP, 1978.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.002}         Max grid cell size
%       tmax        scalar {0.2}           Maximum time
%       tstep       scalar {0.0025}        Time step size
%       ischeme     scalar {2}             Time discretization scheme
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.002;
            'tmax',     0.2;
            'tstep',    0.0025;
            'sfun',     'sflag1';
            'ischeme',  2;
            'iplot',    1;
            'tol',      0.05;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


fea.sdim = { 'x' };
fea.geom.objects = { gobj_line(0,1) };
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );

fea = addphys(fea,@compressibleeuler);

rhoL  = 1;
uL    = 0;
pL    = 1;

rhoR  = 0.125;
uR    = 0;
pR    = 0.1;

x     = fea.grid.p(:);
x0    = (max(x) - min(x))/2;

init0 = { sprintf('%g+(x>%g)*%g',rhoL,x0,rhoR-rhoL), ...
          sprintf('%g+(x>%g)*%g',uL,x0,uR-uL), ...
          sprintf('%g+(x>%g)*%g',pL,x0,pR-pL) };

fea = parsephys(fea);
fea = parseprob(fea);

[fea.sol.u,fea.sol.t] = solvetime( fea, 'init', init0, 'fid', fid, ...
                                   'ischeme', 2, 'tstep', opt.tstep, 'tmax', opt.tmax );


% Error checking.
[r_ref,u_ref, p_ref] = analytic_sod( opt.tmax, x );
r = evalexprp( fea.dvar{1}, fea );
u = evalexprp( fea.dvar{2}, fea );
p = evalexprp( fea.dvar{3}, fea );
out.err = [ sum(abs(r_ref-r))/length(x), ...
            sum(abs(u_ref-u))/length(x), ...
            sum(abs(p_ref-p))/length(x) ];
out.pass = all(out.err<opt.tol);


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(2,2,1)
  plot( x, r, 'k-' )
  hold on
  plot( x, r_ref, 'b--' )
  legend('Computed','Analytic solution')
  title('Density (rho)')
  subplot(2,2,2)
  plot( x, u, 'k-' )
  hold on
  plot( x, u_ref, 'b--' )
  title('Velocity (u)')
  subplot(2,2,3)
  plot( x, p, 'k-' )
  hold on
  plot( x, p_ref, 'b--' )
  title('Pressure (p)')
  subplot(2,2,4)
  plot( x, p./r, 'k-' )
  hold on
  plot( x, p_ref./r_ref, 'b--' )
  title('Non-dimensional temperature (p/rho)')
end


if( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ r, u, p, e ] = analytic_sod( t, x )
% Analytic solution to Sod's Shock Tube problem.

%___|_______|___|_____|_________|_______________
%   x1      x2  x0    x3        x4

% Ref. Gogol https://www.mathworks.com/matlabcentral/fileexchange/46311-sod-shock-tube-problem-solver BSD License.
% Ref. http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html

if( ~nargin )
  t = 0.2;
end
if( nargin<2 )
  x = linspace(0,1,1000);
end
n  = length(x);
x0 = (max(x) - min(x))/2;

r_l = 1;
u_l = 0;
p_l = 1;

u_r = 0;
r_r = 0.125;
p_r = 0.1;

ga  = 1.4;
mu  = sqrt( (ga-1)/(ga+1) );
c_l = sqrt( (ga*p_l/r_l) );
c_r = sqrt( (ga*p_r/r_r) );

sod_func = @(P) (P - p_r)*sqrt((1-mu*mu)*(r_r*(P + mu*mu*p_r))^-1) - ...
           (power(p_l , (ga-1)/(2*ga)) - power(P , (ga-1)/(2*ga))) * ...
           sqrt(((1-mu*mu*mu*mu)*p_l^(1/ga))*(mu*mu*mu*mu*r_l)^-1);

p_post = fzero( sod_func, pi );
r_post = r_r*(( (p_post/p_r) + mu^2 )/(1 + mu*mu*(p_post/p_r)));
u_post = 2*(sqrt(ga)/(ga - 1))*(1 - power(p_post, (ga - 1)/(2*ga)));
u_shock = u_post*((r_post/r_r)/( (r_post/r_r) - 1));
r_middle = (r_l)*power((p_post/p_l),1/ga);

x1 = x0 - c_l*t;
x3 = x0 + u_post*t;
x4 = x0 + u_shock*t;
c_2 = c_l - ((ga - 1)/2)*u_post;
x2 = x0 + (u_post - c_2)*t;

r = zeros(size(x));
u = r;
p = r;

ix = x < x1;
r(ix) = r_l;
u(ix) = u_l;
p(ix) = p_l;

ix = x1 <= x & x <= x2;
c = mu*mu*((x0 - x(ix))/t) + (1 - mu*mu)*c_l;
r(ix) = r_l*power((c/c_l),2/(ga - 1));
u(ix) = (1 - mu*mu)*( (-(x0-x(ix))/t) + c_l);
p(ix) = p_l*power((r(ix)/r_l),ga);

ix = x2 <= x & x <= x3;
r(ix) = r_middle;
u(ix) = u_post;
p(ix) = p_post;

ix = x3 <= x & x <= x4;
r(ix) = r_post;
u(ix) = u_post;
p(ix) = p_post;

ix = x4 < x;
r(ix) = r_r;
u(ix) = u_r;
p(ix) = p_r;

e = p./((ga - 1)*r);
