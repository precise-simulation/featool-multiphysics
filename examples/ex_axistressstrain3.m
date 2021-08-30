function [ fea, out ] = ex_axistressstrain3( varargin )
%EX_AXISTRESSSTRAIN3 Disc with fixed edge and central point load axisymmetric stress-strain.
%
%   [ FEA, OUT ] = EX_AXISTRESSSTRAIN3( VARARGIN ) Example to calculate displacements and stresses
%   in a disc with fixed outer edge with a central point load in axisymmetric/cylindrical coordinates.
%
%   Ref. Chapter 7. Benchmark 3: Point Loaded SS Circular Plate Bending.
%   [1] Advanced Finite Element Methods (ASEN 6367) - Spring 2013. Department of Aerospace
%       Engineering Sciences University of Colorado at Boulder.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       p           scalar {1e3}           Load force
%       E           scalar {1e3}           Modulus of elasticity
%       nu          scalar {1/3}           Poissons ratio
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'p',        1e3;
            'E',        1e3;
            'nu',       1/3;
            'sfun',     'sflag1';
            'iplot',    1;
            'igrid',    1;
            'tol',      0.02;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid.
a = 0;
b = 10;
fea.grid = rectgrid( 20, 20, [a b;-0.5 0.5] );
if( opt.igrid<0 )
  fea.grid = quad2tri( fea.grid );
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
fea.phys.css.bdr.coef{1,5} = bctype;

bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef{2,4} = -opt.p/(2*pi);
fea.phys.css.bdr.coef{1,end} = bccoef;


% Set point constraint w=0 at (b,0).
[~,ix] = min( (fea.grid.p(1,:)-b).^2 + (fea.grid.p(2,:)-0).^2 );
fea.pnt(1).index = ix;
fea.pnt(1).type  = 'constr';
fea.pnt(1).dvar  = 2;
fea.pnt(1).expr  = 0;


% Solve.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')), 'fid', fid );


% Postprocessing.
D = opt.E*1^3/(12*(1-opt.nu^2));
u_ref_r = [num2str(-opt.p),'/(8*pi*',num2str(D),')*(',num2str((3+opt.nu)/(1+opt.nu)-1),'-2*log(r/',num2str(b),'))*r*z'];
u_ref_z = [num2str(-opt.p),'/(16*pi*',num2str(D),')*(',num2str((3+opt.nu)/(1+opt.nu)),'*(',num2str(b^2),'-r^2)+2*r^2*log(r/',num2str(b),'))'];
if( opt.iplot>0 )
  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'r*u' )
  title('computed r-displacement')

  subplot(2,2,2)
  postplot( fea, 'surfexpr', u_ref_r )
  title('exact r-displacement')

  subplot(2,2,3)
  postplot( fea, 'surfexpr', 'w' )
  title('computed z-displacement')

  subplot(2,2,4)
  postplot( fea, 'surfexpr', u_ref_z )
  title('exact z-displacement')
end


% Error checking.
u_r = evalexpr( 'r*u', fea.grid.p, fea );
u_z = evalexpr( 'w', fea.grid.p, fea );
u_ref_r = evalexpr( u_ref_r, fea.grid.p, fea );
u_ref_z = evalexpr( u_ref_z, fea.grid.p, fea );
ix = find( ~isnan(u_ref_r + u_ref_z) );
out.err(1) = norm( u_ref_r(ix) - u_r(ix) )/norm( u_ref_r(ix) );
out.err(2) = norm( u_ref_z(ix) - u_z(ix) )/norm( u_ref_z(ix) );
out.pass = all( out.err < opt.tol );


if( nargout==0 )
  clear fea out
end
