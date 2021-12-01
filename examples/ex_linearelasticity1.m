function [ fea, out ] = ex_linearelasticity1( varargin )
%EX_LINEARELASTICITY1 Example for solid stress-strain on a cube.
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY1( VARARGIN ) Example to
%   calculate displacements and stresses for a cube stretched in
%   axis-aligned directions.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {1}             Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       force       scalar {1}             Load force
%       idir        scalar 1/{2,3}         Stress direction
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Max grid cell size
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'E',        1;
            'nu',       0.3;
            'force',    1;
            'idir',     1;
            'igrid',    0;
            'hmax',     0.1;
            'sfun',     'sflag1';
            'iplot',    1;
            'psecheck', 0;
            'tol',      sqrt(eps);
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.geom.objects = { gobj_block() };
fea.sdim = { 'x' 'y' 'z' };   % Coordinate names.


% Grid generation.
switch opt.igrid
  case -1
    fea.grid = blockgrid( round(1/opt.hmax) );
    fea.grid = hex2tet( fea.grid );
  case 0
    fea.grid = blockgrid( round(1/opt.hmax) );
  case 1
    fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Boundary conditions.
s = sign(opt.idir); opt.idir = abs(opt.idir);
dtol     = opt.hmax;
fixbdr   = findbdr( fea, [fea.sdim{opt.idir},'<',num2str(dtol)] );     % "Right" boundary number.
forcebdr = findbdr( fea, [fea.sdim{opt.idir},'>',num2str(1-dtol)] );   % "Left" boundary number.


% Problem definition.
fea = addphys(fea,@linearelasticity);
fea.phys.el.eqn.coef{1,end} = { opt.nu };
fea.phys.el.eqn.coef{2,end} = { opt.E  };
fea.phys.el.sfun            = { opt.sfun opt.sfun opt.sfun };

% Fix first boundary (set idir displacement to zero, Dirichlet BC).
bctype = mat2cell( zeros(3,n_bdr), [1 1 1], ones(1,n_bdr) );
bctype{opt.idir,fixbdr} = 1;
fea.phys.el.bdr.coef{1,5} = bctype;

% Apply load to second opposite boundary.
bccoef = mat2cell( zeros(3,n_bdr), [1 1 1], ones(1,n_bdr) );
bccoef{opt.idir,forcebdr} = s*opt.force;
fea.phys.el.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', fid, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')) );


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', fea.dvar{opt.idir} )
  title( [fea.sdim{opt.idir},'-displacement'] )
end


% Error checking.
x = linspace(0,1,7);
x = x(2:end-1)';
[x,y,z] = ndgrid(x,x,x);
p = [x(:) y(:) z(:)]';
ui = evalexpr( fea.dvar{opt.idir}, p, fea );
ui_ref = s*p(opt.idir,:)';
out.ui_err = norm( ui - ui_ref );
out.pserr  = [];
out.peerr  = [];
for i=1:6
  out.serr(i) = norm( abs( evalexpr( fea.phys.el.eqn.vars{5+i, end}, p, fea ) - s*(i==opt.idir) ) );
  out.eerr(i) = norm( abs( evalexpr( fea.phys.el.eqn.vars{14+i,end}, p, fea ) - s*(i==opt.idir) + s*0.3*(i<=3 & i~=opt.idir)) ) ;
  if( opt.psecheck && i<=3 )
    out.pserr(i) = norm( abs( evalexpr( fea.phys.el.eqn.vars{11+i,end}, p, fea ) - (i==1) ) );
    out.peerr(i) = norm( abs( evalexpr( fea.phys.el.eqn.vars{20+i,end}, p, fea ) - (i==1) + 0.3*(i~=1) ) );
  end
end
out.serr = [ out.serr norm( abs( evalexpr( fea.phys.el.eqn.vars{1, end}, p, fea ) - 1 ) ) ];
out.pass = all( [out.ui_err out.serr out.eerr out.pserr out.peerr] < opt.tol );


if ( nargout==0 )
  clear fea out
end
