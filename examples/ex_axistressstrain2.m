function [ fea, out ] = ex_axistressstrain2( varargin )
%EX_AXISTRESSSTRAIN2 Example for a pressurized hollow sphere axisymmetric stress-strain.
%
%   [ FEA, OUT ] = EX_AXISTRESSSTRAIN2( VARARGIN ) Example to calculate displacements and stresses
%   in a pressurized hollow sphere in axisymmetric/cylindrical coordinates.
%
%   Ref. 4.1.4 Pressurized hollow sphere.
%   [1] Applied Mechanics of Solids, Allan F. Bower, 2012 (http://solidmechanics.org/).
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       a           scalar {1}             Cylinder inner radius
%       b           scalar {2}             Cylinder outer radius
%       p           scalar {20e4}          Load force
%       E           scalar {200e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       igrid       scalar 0/{<0}          Cell type (0=quadrilaterals, <0=triangles)
%       sfun        string {sflag2}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'a',        1;
            'b',        2;
            'p',        20e4;
            'E',        200e9;
            'nu',       0.3;
            'igrid',    0;
            'sfun',     'sflag2';
            'iplot',    1;
            'tol',      5e-3;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid.
a = opt.a;
b = opt.b;
if ( opt.igrid==1 )
  error('ex_axistressstrain2: unstructured grid not supported.')
else
  fea.grid = ringgrid( 12, 72, a, b );
  fea.grid = delcells( fea.grid, selcells( fea.grid, '(x<=eps) | (y<=eps)') );
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
bctype{1,4} = 1;
bctype{2,3} = 1;
fea.phys.css.bdr.coef{1,5} = bctype;

bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef{1,1} = ['-r*nr*',num2str(opt.p),];
bccoef{2,1} = ['-r*nz*',num2str(opt.p),];
fea.phys.css.bdr.coef{1,end} = bccoef;


% Solve.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')), 'fid', fid );


% Postprocessing.
n = 20;
r = linspace(a,b,n);
z = zeros(1,n);
u_ref = 1./(2*opt.E*(b^3-a^3)*r'.^2) .* (2*(opt.p*a^3)*(1-2*opt.nu)*r'.^3+opt.p*(1+opt.nu)*b^3*a^3);
u = evalexpr( 'r*u', [r;z], fea );
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt((r*u)^2+w^2)', 'arrowexpr', {'r*u' 'w'} )
  title('computed displacement')
  subplot(1,2,2), hold on
  plot(u_ref,r,'r-')
  plot(u,r,'b.')
  title( 'radial displacement')
  legend('exact solution','computed solution')
  xlabel('r')
  grid on
end


% Error checking.
out.err  = norm( u_ref - u )/norm( u_ref );
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
