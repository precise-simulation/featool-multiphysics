function [ fea, out ] = ex_planestress5( varargin )
%EX_PLANESTRESS5 Plane stress example for an elliptic membrane.
%
%   [ FEA, OUT ] = EX_PLANESTRESS5( VARARGIN ) Example to calculate displacements and stresses
%   for an elliptic membrane with a hole in it. NAFEMS benchmark example LE1.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {210e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       hmax        scalar {0.1}           Max grid cell size
%       sfun        string {sflag2}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        210e9; ...
  'nu',       0.3; ...
  'hmax',     0.1; ...
  'sfun',     'sflag2'; ...
  'iplot',    1; ...
  'igeom',    1; ...
  'tol',      0.05; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
gobj1 = gobj_ellipse( [0 0], 3.25, 2.75, 'E1' );
gobj2 = gobj_ellipse( [0 0], 2, 1, 'E2' );
gobj3 = gobj_rectangle( -3.25, 3.25, -2.75, 0, 'R1' );
gobj4 = gobj_rectangle( -3.25, 0, 0, 2.75, 'R2' );
fea.geom.objects = { gobj1 gobj2 gobj3 gobj4 };
if( opt.igeom==1 )
  fea = geom_apply_formula( fea, 'E1-E2-R1-R2' );
else
  fea = geom_apply_formula( fea, 'E1-E2' );
  fea = geom_apply_formula( fea, 'CS1-R1' );
  fea = geom_apply_formula( fea, 'CS2-R2' );
end
fea.sdim = { 'x' 'y' };


% Grid generation.
if( opt.igeom==1 )
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid, 'gridgen', 'gridgen2d' );
else
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Boundary conditions.
dtol  = sqrt(eps);
lbdr  = findbdr( fea, ['x<=',num2str(dtol)] );   % Left boundary number.
lobdr = findbdr( fea, ['y<=',num2str(dtol)] );   % Lower boundary number.


% Problem definition.
E11 = opt.E/(1-opt.nu^2);
E12 = opt.nu*E11;
E22 = E11;
E33 = opt.E/(1+opt.nu)/2;

fea = addphys(fea,@planestress);      % Add plane stress physics mode.
fea.phys.pss.eqn.coef{1,end} = { opt.nu };
fea.phys.pss.eqn.coef{2,end} = { opt.E  };
fea.phys.pss.sfun            = { opt.sfun opt.sfun };   % Set shape functions.

bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bctype{1,lbdr}  = 1;
bctype{2,lobdr} = 1;
fea.phys.pss.bdr.coef{1,5} = bctype;

% Add normal load to outer boundary.
dtol = 1e-3;
i_o = findbdr( fea, ['sqrt(x.^2+y.^2)>=',num2str(2.75-dtol)] );
bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef{1,i_o} = 'nx*10e6';
bccoef{2,i_o} = 'ny*10e6';
fea.phys.pss.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea       = parsephys(fea);
fea       = parseprob(fea);
fea.sol.u = solvestat( fea, 'fid', fid );


% Postprocessing.
s_sy = [num2str(E12),'*ux+',num2str(E11),'*vy'];
if ( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_sy, 'isoexpr', s_sy )
  title('Stress, x-component')
end


% Error checking.
sy_D     = evalexpr( s_sy, [2;0]+sqrt(eps)*1e1, fea );
out.sy_D = sy_D;
out.err  = abs(sy_D - 92.7e6)/92.7e6;
out.pass = out.err <= opt.tol;


if ( nargout==0 )
  clear fea out
end
