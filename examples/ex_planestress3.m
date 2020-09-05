function [ fea, out ] = ex_planestress3( varargin )
%EX_PLANESTRESS3 NAFEMS benchmark model for stress of a tapered membrane.
%
%   [ FEA, OUT ] = EX_PLANESTRESS3( VARARGIN ) NAFEMS benchmark for stress of a tapered membrane.
%   Two test cases are modeled, the first with the x displacements fixes on the left boundary while
%   a horizontal force is applied to the right boundary. Secondly a gravitational body force is applied
%   while the left boundary is fixed. Reference: Linear Statics Benchmarks Vol. 1, NAFEMS Ltd., 1987.
%
%   Accepts the following property/value pairs.
%
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {210e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       thick       scalar {0.1}           Plate thickness
%       hmax        scalar {1/20}          Max grid cell size
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2020 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        210e9; ...
  'nu',       0.3; ...
  'thick',    0.1; ...
  'icase',    1; ...
  'hmax',     0.3; ...
  'sfun',     'sflag2'; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj1 = gobj_polygon( [0 4 4 0 0;0 1 3 4 0]', 'P1' );
fea.geom.objects = { gobj1 };
fea.sdim = { 'x' 'y' };


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );


% Add plane stress physics mode.
fea = addphys(fea,@planestress);
fea.phys.pss.eqn.coef{1,end} = { opt.nu };
fea.phys.pss.eqn.coef{2,end} = { opt.E  };
fea.phys.pss.sfun            = { opt.sfun opt.sfun };
if( opt.icase==2 )
  fea.phys.pss.eqn.coef{4,end} = { -9.81*7000 };
end


% Set boundary conditions.
dtol = 0.1;
lbdr = findbdr( fea, ['x<',num2str(dtol)] );     % Left boundary number.
rbdr = findbdr( fea, ['x>',num2str(1-dtol)] );   % Right boundary number.
n_bdr  = max(fea.grid.b(3,:));        % Number of boundaries.
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
if( opt.icase==1 )
  bctype{1,lbdr} = 1;
  bccoef{1,rbdr} = 1e7/opt.thick;
else
  bctype{1,lbdr} = 1;
  bctype{2,lbdr} = 1;
end
fea.phys.pss.bdr.coef{1,end} = bccoef;
fea.phys.pss.bdr.coef{1,5}   = bctype;


% Parse and solve problem.
fea       = parsephys(fea);             % Check and parse physics modes.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
if( opt.icase==1 )
  s_title = fea.phys.pss.eqn.vars{5,1};
  s_expr  = fea.phys.pss.eqn.vars{5,2};
else
  s_title = fea.phys.pss.eqn.vars{7,1};
  s_expr  = fea.phys.pss.eqn.vars{7,2};
end
if ( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_expr, 'isoexpr', s_expr )
  title( s_title )
end


% Error checking.
s_02  = evalexpr( s_expr, [0;2], fea );
s_ref = (61.3e6*(opt.icase==1)) + (-0.2e6*(opt.icase==2));
out.stress = s_02;
out.err    = abs(s_02-s_ref)/s_ref;
out.pass   = out.err<0.01;


if ( nargout==0 )
  clear fea out
end
