function [ fea, out ] = ex_planestress1( varargin )
%EX_PLANESTRESS1 Example for plane stress for hole in plate.
%
%   [ FEA, OUT ] = EX_PLANESTRESS1( VARARGIN ) Example to calculate displacements and stresses
%   for a hole in plate configuration under plane stress assumption.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {210e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       diam        scalar {0.01}          Diameter of hole
%       thick       scalar {0.001}         Plate thickness
%       force       scalar {1000}          Load force
%       sx_ref      scalar {3e7}           Reference stress in x-direction
%       hmax        scalar {0.00125}       Max grid cell size
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
  'E',        210e9; ...
  'nu',       0.3; ...
  'diam',     0.01; ...
  'thick',    0.001; ...
  'force',    1000; ...
  'sx_ref',   3e7; ...
  'hmax',     0.00125; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'tol',      0.05; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model, geometry, and grid parameters.
h         = 0.05;            % Height of 1/4 rectangular domain.
l         = 0.05;            % Length of 1/4 rectangular domain.
xc        = 0;               % x-coordinate of hole center.
yc        = 0;               % y-coordinate of hole center.
area      = 2*h*opt.thick;   % Area on which load force is applied.


% Geometry definition.
gobj1 = gobj_rectangle( 0, l, 0, h, 'R1' );
gobj2 = gobj_circle( [xc,yc], opt.diam/2, 'C1' );
fea.geom.objects = { gobj1 gobj2 };
fea = geom_apply_formula( fea, 'R1-C1' );
fea.sdim = { 'x' 'y' };   % Coordinate names.


% Grid generation.
fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
n_bdr    = max(fea.grid.b(3,:));        % Number of boundaries.


% Boundary conditions.
dtol  = opt.diam/10;
lbdr  = findbdr( fea, ['x<=',num2str(dtol)] );     % Left boundary number.
rbdr  = findbdr( fea, ['x>=',num2str(l-dtol)] );   % Right boundary number.
lobdr = findbdr( fea, ['y<=',num2str(dtol)] );     % Lower boundary number.


% Problem definition.
E11 = opt.E/(1-opt.nu^2);
E12 = opt.nu*E11;
E22 = E11;
E33 = opt.E/(1+opt.nu)/2;

if ( opt.iphys==1 )

  fea = addphys(fea,@planestress);      % Add plane stress physics mode.
  fea.phys.pss.eqn.coef{1,end} = { opt.nu };
  fea.phys.pss.eqn.coef{2,end} = { opt.E  };
  fea.phys.pss.sfun            = { opt.sfun opt.sfun };   % Set shape functions.

  bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
  bctype{1,lbdr}  = 1;
  bctype{2,lobdr} = 1;
  fea.phys.pss.bdr.coef{1,5}   = bctype;
  bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
  bccoef{1,rbdr} = opt.force/area;
  fea.phys.pss.bdr.coef{1,end} = bccoef;

  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u' 'v' };                  % Dependent variable names.
  fea.sfun  = { opt.sfun opt.sfun };        % Shape functions.

  % Define equation system.
  fea.eqn.a.form = { [2 3;2 3] [3 2;2 3];   ...
                     [3 2;2 3] [2 3;2 3] }; ...
  fea.eqn.a.coef = { {E11 E33} {E12 E33};   ...
                     {E33 E12} {E33 E11} };
  fea.eqn.f.form = { 1 1 };
  fea.eqn.f.coef = { 0 0 };

  % Define boundary conditions.
  fea.bdr.d     = cell(2,n_bdr);
  fea.bdr.n     = cell(2,n_bdr);
 [fea.bdr.n{:}] = deal(0);              % Assign zero to all Neumann boundaries.

  fea.bdr.n{1,rbdr}  = opt.force/area;  % Set horizontal load force on right boundary.

  fea.bdr.d{1,lbdr}  = 0;               % Set zero horizontal displacement on left boundary.

  fea.bdr.d{2,lobdr} = 0;               % Set zero vertical displacement on lower boundary.

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
s_sx = [num2str(E11),'*ux+',num2str(E12),'*vy'];
if ( opt.iplot>0 )
  figure
  postplot(fea,'surfexpr',s_sx,'isoexpr',s_sx,'isomap','k')
  title('Stress, x-component')
end


% Error checking.
[~,sx_max] = minmaxsubd( s_sx, fea );
out.sx_max = sx_max;
out.err    = opt.sx_ref-sx_max;
out.pass   = (sx_max>(opt.sx_ref*(1-opt.tol)))&&(sx_max<(opt.sx_ref*(1+opt.tol)));


if ( nargout==0 )
  clear fea out
end
