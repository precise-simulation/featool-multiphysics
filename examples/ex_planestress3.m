function [ fea, out ] = ex_planestress3( varargin )
%EX_PLANESTRESS3 NAFEMS benchmarks IC1-4 linear static stress analysis of a tapered membrane.
%
%   [ FEA, OUT ] = EX_PLANESTRESS3( VARARGIN ) NAFEMS benchmarks IC1-4
%   for linear static plane stress analysis of a tapered membrane.
%   Four test cases are modeled, the first with a horizonal load on
%   the left edge, second with horizonal volume force, third with a
%   vertical shear load on the left edge, and fourth with a vertical
%   volume (gravity) force.
%
%   Reference: Linear Statics Benchmarks Vol. 1, NAFEMS Ltd., 1987.
%
%   Accepts the following property/value pairs.
%
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       icase       scalar 1-4/{1}
%       hmax        scalar {0.3}           Max grid cell size
%       sfun        string {sflag2}        Shape function for displacements
%       solver      string {}              Solver selection default, fenics
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'E',     210e9;
            'nu',    0.3;
            'thick', 0.1;
            'icase', 1;
            'hmax',  0.3;
            'sfun',  'sflag2';
            'solver',   '';
            'iplot', 1;
            'tol',   0.01;
            'fid',   1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_polygon( [0 4 4 0 0 0;0 1 3 4 2 0]', 'P1' );
fea.geom.objects = { gobj };
fea.sdim = { 'x' 'y' };


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );


% Add plane stress physics mode.
fea = addphys(fea,@planestress);
fea.phys.pss.eqn.coef{1,end} = { opt.nu };
fea.phys.pss.eqn.coef{2,end} = { opt.E  };
fea.phys.pss.sfun            = { opt.sfun opt.sfun };
if( opt.icase == 2 )
  fea.phys.pss.eqn.coef{4,end} = { 9.81*7000 };
elseif( opt.icase == 4 )
  fea.phys.pss.eqn.coef{5,end} = { -9.81*7000 };
end


% Set boundary conditions.
dtol = 0.1;
lbdr = findbdr( fea, ['x<',num2str(dtol)] );     % Left boundary number.
rbdr = findbdr( fea, ['x>',num2str(1-dtol)] );   % Right boundary number.
n_bdr  = max(fea.grid.b(3,:));        % Number of boundaries.
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
switch( opt.icase )
  case {1,2}
    [bctype{1,lbdr}] = deal(1);
    if( opt.icase == 1 )
      bccoef{1,rbdr} = 1e7/opt.thick;
    end
    fea.pnt(1).index = [0;2];
    fea.pnt(1).type  = 'constr';
    fea.pnt(1).dvar  = 'u';
    fea.pnt(1).expr  = 0';
    fea.pnt(2).index = [0;2];
    fea.pnt(2).type  = 'constr';
    fea.pnt(2).dvar  = 'v';
    fea.pnt(2).expr  = 0';
  case 3
    [bctype{:,lbdr}] = deal(1);
    bccoef{2,rbdr} = 1e7/opt.thick;
  case 4
    [bctype{:,lbdr}] = deal(1);
end
fea.phys.pss.bdr.coef{1,end} = bccoef;
fea.phys.pss.bdr.coef{1,5}   = bctype;


% Parse and solve problem.
fea       = parsephys(fea);             % Check and parse physics modes.
fea       = parseprob(fea);             % Check and parse problem struct.
if( strcmp(opt.solver,'fenics') )
  fea = fenics(fea,'fid',fid);
else
  fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.
end


% Postprocessing.
if( opt.icase <= 2 )
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
s_02 = evalexpr( s_expr, [0;2], fea );
s_ref = [61.3, 0.247, 26.9, -0.2]*1e6;
out.stress = s_02;
out.err    = abs(s_02 - s_ref(opt.icase))/s_ref(opt.icase);
out.pass   = out.err < opt.tol;


if ( nargout==0 )
  clear fea out
end
