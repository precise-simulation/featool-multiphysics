function [ fea, out ] = ex_planestress4( varargin )
%EX_PLANESTRESS4 Example of thermally induced stress.
%
%   [ FEA, OUT ] = EX_PLANESTRESS4( VARARGIN ) NAFEMS T1 Benchmark
%   example for thermally induced stress. A 20x20 mm plate is
%   subjected to thermal load at a circular spot with radius 1
%   mm. Plane stress and symmetry can be assumed. The stress in the
%   y-direction is sought just outside the thermal spot.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {100e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       hmax        scalar {4e-4}          Max grid cell size
%       sfun        string {sflag1}        Shape function for displacements
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        100e9; ...
  'nu',       0.3; ...
  'sy_ref',   50e6; ...
  'hmax',     3e-4; ...
  'sfun',     'sflag2'; ...
  'igeom',    1; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'tol',      0.1; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Model, geometry, and grid parameters.
l  = 1e-2;   % Length of quarter domain.
r  = 1e-3;   % Radius of thermal spot.
xc = 0;      % x-coordinate of spot center.
yc = 0;      % y-coordinate of spot center.


% Geometry definition.
switch( opt.igeom )
  case 1
    gobj1 = gobj_rectangle( 0, l, 0, l, 'R1' );
    gobj2 = gobj_circle( [0 0], r, 'C1' );
    gobj3 = gobj_rectangle( 0, l, 0, l, 'R2' );
    gobj4 = gobj_circle( [0 0], r, 'C2' );
    fea.geom.objects = { gobj1 gobj2 gobj3 gobj4 };
    fea.geom = geom_apply_formula( fea.geom, 'R1-C1' );
    fea.geom = geom_apply_formula( fea.geom, 'R2&C2' );
  case 2
    gobj1 = gobj_rectangle( 0, l, 0, l, 'R1' );
    gobj2 = gobj_circle( [0 0], r, 'C1' );
    gobj3 = gobj_rectangle( 0, l, 0, l, 'R2' );
    fea.geom.objects = { gobj1 gobj2 gobj3 };
    fea.geom = geom_apply_formula( fea.geom, 'R1&C1' );
  case 3
    gobj1 = gobj_rectangle( 0, l, 0, l );
    gobj2 = gobj_circle(  [0 0], r );
    gobj3 = gobj_polygon( [0 0; r 0;r -r; -r -r; -r r; 0 r] );
    fea.geom.objects = { gobj1 gobj2 gobj3 };
    fea.geom = geom_apply_formula( fea.geom, 'C1-P1' );
end
fea.sdim = { 'x' 'y' };   % Coordinate names.


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
n_bdr    = max(fea.grid.b(3,:));   % Number of boundaries.


% Boundary conditions.
dtol  = r/5;
lbdr  = findbdr( fea, ['x<=',num2str(dtol)] );   % Left boundary number.
lobdr = findbdr( fea, ['y<=',num2str(dtol)] );   % Lower boundary number.


% Problem definition.
if( opt.iphys==1 )

  fea = addphys( fea, @planestress );
  fea.phys.pss.eqn.coef{1,end} = { opt.nu };
  fea.phys.pss.eqn.coef{2,end} = { opt.E  };
  if( sum(fea.grid.s==1) < sum(fea.grid.s==2) )
    fea.phys.pss.eqn.coef{6,end} = { 1e-3, 0 };
  else
    fea.phys.pss.eqn.coef{6,end} = { 0, 1e-3 };
  end
  fea.phys.pss.eqn.coef{7,end} = { 1 };
  fea.phys.pss.sfun            = { opt.sfun opt.sfun };

  bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
  bctype{1,lbdr(1)}  = 1;
  bctype{1,lbdr(2)}  = 1;
  bctype{2,lobdr(1)} = 1;
  bctype{2,lobdr(2)} = 1;
  fea.phys.pss.bdr.coef{1,5} = bctype;

  fea = parsephys(fea);

  s_sy = fea.phys.pss.eqn.vars{6,end};

else

  E11 = opt.E/(1-opt.nu^2);
  E12 = opt.nu*E11;
  E22 = E11;
  E33 = opt.E/(1+opt.nu)/2;
  if( sum(fea.grid.s==1) < sum(fea.grid.s==2) )
    fea.expr = { 'alfaT', { opt.E/(1-opt.nu)*1e-3 0 } };
  else
    fea.expr = { 'alfaT', { 0 opt.E/(1-opt.nu)*1e-3 } };
  end

  fea.dvar  = { 'u' 'v' };
  fea.sfun  = { opt.sfun opt.sfun };

  % Define equation system.
  fea.eqn.a.form = { [2 3;2 3] [3 2;2 3];   ...
                     [3 2;2 3] [2 3;2 3] }; ...
  fea.eqn.a.coef = { {E11 E33} {E12 E33};   ...
                     {E33 E12} {E33 E11} };
  fea.eqn.f.form = { 2 3 };
  fea.eqn.f.coef = { 'alfaT' 'alfaT' };

  % Define boundary conditions.
  fea.bdr.d     = cell(2,n_bdr);
  fea.bdr.n     = cell(2,n_bdr);
 [fea.bdr.n{:}] = deal(0);

  fea.bdr.d{1,lbdr(1)}  = 0;
  fea.bdr.d{1,lbdr(2)}  = 0;
  fea.bdr.d{2,lobdr(1)} = 0;
  fea.bdr.d{2,lobdr(2)} = 0;

  s_sy = [num2str(E12),'*ux + ',num2str(E11),'*vy - alfaT'];
end

% Parse and solve problem.
fea       = parseprob( fea );               % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'fid', fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_sy, 'isoexpr', s_sy, 'isolev', 50 )
  title( 'Stress, y-component' )
end


% Error checking.
sy_D     = evalexpr( s_sy, [1e-3;1e-6], fea );
out.sy_D = sy_D;
out.err  = abs(sy_D - opt.sy_ref)/opt.sy_ref;
out.pass = out.err <= opt.tol;


if( nargout==0 )
  clear fea out
end
