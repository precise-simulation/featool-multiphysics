function [ fea, out ] = ex_convdiff5( varargin )
%EX_CONVDIFF5 2D Convection and diffusion equation with high Peclet number.
%
%   [ FEA, OUT ] = EX_CONVDIFF5( VARARGIN ) Convection and diffusion equation on
%   a unit square with high Peclet (Cell Reynolds) number requiring artificial
%   stabilization/numerical diffusion. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar {0}/1           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {1/20}          Max grid cell size
%       a           scalar {cos(pi/3)}     Convection velocity in x-direction
%       b           scalar {sin(pi/3)}     Convection velocity in y-direction
%       cd          scalar {1e-4}          Diffusion coefficient
%       sfun        string {sflag1}        Shape function
%       artstab     scalar {1}/0/2         Artificial stabilization (0=no stabilization)
%                                          1=isotropic diffusion, 2=streamline diffusion
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'igrid',    0; ...
  'hmax',     1/20; ...
  'a',        cos(pi/3); ...
  'b',        sin(pi/3); ...
  'cd',       1e-4; ...
  'sfun',     'sflag1'; ...
  'artstab',  1; ...
  'iplot',    1; ...
  'tol',      0.1; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition and grid generation.
fea.geom.objects = { gobj_rectangle() };
switch opt.igrid
  case -1
    fea.grid = rectgrid( round(1/opt.hmax) );
    fea.grid = quad2tri( fea.grid);
  case 0
    fea.grid = rectgrid( round(1/opt.hmax) );
  case 1
    fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid, 'dprim', false );
end


% Problem definition.
fea.sdim = { 'x' 'y' };                     % Coordinate names.
fea = addphys( fea, @convectiondiffusion ); % Add convection and diffusion physics mode.
fea.phys.cd.sfun            = { opt.sfun }; % Set shape function.
fea.phys.cd.eqn.coef{2,4}   = { opt.cd };   % Set diffusion coefficient.
fea.phys.cd.eqn.coef{3,4}   = { opt.a  };   % Convection velocity in x-direction.
fea.phys.cd.eqn.coef{4,4}   = { opt.b  };   % Convection velocity in y-direction.
fea.phys.cd.eqn.coef{5,4}   = { 1 };
fea.phys.cd.bdr.sel         = [1 1 1 1];
fea.phys.cd.bdr.coef{1,end} = {1 0 0 1};


% Numerical stabilization.
switch opt.artstab
  case 1
    fea.phys.cd.prop.artstab.id      = 1;
    fea.phys.cd.prop.artstab.id_coef = 0.5;
  case 2
    fea.phys.cd.prop.artstab.sd      = 1;
    fea.phys.cd.prop.artstab.sd_coef = 0.25;
end


% Parse and solve problem.
fea       = parsephys( fea );               % Check and parse physics modes.
fea       = parseprob( fea );               % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'fid', fid );   % Call to stationary solver.


% Postprocessing.
x = linspace( 0, 1, 2*(1/opt.hmax)+1 );
y = 0.8*ones(size(x));
c0p8 = evalexpr( 'c', [x;y], fea );
c_ref = 1 + 0.1/0.05*x';
c_ref(c_ref>1.9275) = 1.9275;
if( opt.iplot>0 )
  figure
  subplot( 1, 2, 1 )
  postplot( fea, 'surfexpr', 'c', 'isoexpr', 'c' )
  title( 'Solution c' )
  subplot( 1, 2, 2 )
  plot( x, c_ref, 'r' )
  hold on
  plot( x, c0p8, 'b--' )
  title( 'Solution at y = 0.8' )
  xlabel( 'x' )
  ylabel( 'c' )
  grid on
end


% Error checking.
ix = find(x<=0.9);
out.err  = norm( c_ref(ix) - c0p8(ix) )./norm( c_ref(ix) );
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
