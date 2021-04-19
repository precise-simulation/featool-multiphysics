function [ fea, out ] = ex_axistressstrain1( varargin )
%EX_AXISTRESSSTRAIN1 Example for hollow cylider axisymmetric stress-strain.
%
%   [ FEA, OUT ] = EX_AXISTRESSSTRAIN1( VARARGIN ) Example to calculate displacements and stresses
%   in a hollow cylinder in axisymmetric/cylindrical coordinates.
%
%   Ref. 4.1.9 Long (generalized plane strain) cylinder subjected to internal and external pressure.
%   [1] Applied Mechanics of Solids, Allan F. Bower, 2012 (http://solidmechanics.org/).
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       a           scalar {1.5}           Cylinder inner radius
%       b           scalar {2}             Cylinder outer radius
%       l           scalar {3}             Cylinder length
%       pa          scalar {5e3}           Inner load force
%       pb          scalar {20e4}          Outer load force
%       E           scalar {200e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Max grid cell size
%       sfun        string {sflag2}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'a',        0.9;
            'b',        2;
            'l',        3;
            'pa',       5e3;
            'pb',       20e4;
            'E',        200e9;
            'nu',       0.3;
            'igrid',    0;
            'hmax',     0.01;
            'sfun',     'sflag2';
            'iplot',    1;
            'tol',      5e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid.
a = opt.a;
b = opt.b;
l = opt.l;
fea.geom.objects = { gobj_rectangle( a, b, 0, l, 'R1' ) };
if ( opt.igrid==1 )
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
else
  fea.grid = rectgrid( ceil((b-a)/opt.hmax), ceil(l/opt.hmax), [a b;0 l] );
  if( opt.igrid<0 )
    fea.grid = quad2tri( fea.grid );
  end
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Axisymmetric stress-strain equation definitions.
fea.sdim = { 'r', 'z' };
fea = addphys( fea, @axistressstrain );
fea.phys.css.eqn.coef{1,end} = { opt.nu };
fea.phys.css.eqn.coef{2,end} = { opt.E  };
fea.phys.css.sfun            = { opt.sfun opt.sfun };   % Set shape functions.

% Boundary conditions.
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
[bctype{2,:}] = deal( 1 );
fea.phys.css.bdr.coef{1,5} = bctype;

bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef{1,2} = opt.pb*b;
bccoef{1,4} = opt.pa*a;
fea.phys.css.bdr.coef{1,end} = bccoef;


% Solve.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')), 'fid', fid );


% Postprocessing.
n = 20;
r = linspace(a,b,n);
z = l/2*ones(1,n);
pa = opt.pa; pb = -opt.pb; E = opt.E; nu = opt.nu;
% 4.1.9 http://solidmechanics.org/Text/Chapter4_1/Chapter4_1.php#Sect4_1_9
u_ref  = (1+nu)*a^2*b^2/E/(b^2-a^2)*( (pa-pb)./r' + (1-2*nu)*(pa*a^2 - pb*b^2)/a^2/b^2*r' );
sr_ref = (pa*a^2-pb*b^2)/(b^2-a^2) - a^2*b^2/(b^2-a^2)./(r').^2*(pa-pb);
st_ref = (pa*a^2-pb*b^2)/(b^2-a^2) + a^2*b^2/(b^2-a^2)./(r').^2*(pa-pb);
sz_ref = 2*nu*(pa*a^2-pb*b^2)/(b^2-a^2);
u  = evalexpr( fea.phys.css.eqn.vars{3,2}, [r;z], fea );
w  = evalexpr( fea.phys.css.eqn.vars{4,2}, [r;z], fea );
sr = evalexpr( fea.phys.css.eqn.vars{5,2}, [r;z], fea );
st = evalexpr( fea.phys.css.eqn.vars{6,2}, [r;z], fea );
sz = evalexpr( fea.phys.css.eqn.vars{7,2}, [r;z], fea );
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', fea.phys.css.eqn.vars{3,2} )
  title('r-displacement')
  subplot(1,2,2), hold on
  plot(u_ref,r,'r-')
  plot(u,r,'b.')
  legend('exact solution','computed solution')
  xlabel('r')
  grid on
end


% Error checking.
out.erru  = norm( u_ref - u )/norm( u_ref );
out.errw  = norm( w );
out.errsr = norm( sr_ref - sr )/norm( sr_ref );
out.errst = norm( st_ref - st )/norm( st_ref );
out.errsz = norm( sz_ref - sz )/norm( sz_ref );
out.err   = [ out.erru, out.errw, out.errsr, out.errst, out.errsz ];
out.pass  = all(out.err < opt.tol);


if( nargout==0 )
  clear fea out
end
