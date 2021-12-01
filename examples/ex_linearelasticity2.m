function [ fea, out ] = ex_linearelasticity2( varargin )
%EX_LINEARELASTICITY2 Example for deflection of a bracket.
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY2( VARARGIN ) Example to calculate displacements and
%   stresses for a bracket with a circular hole.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {200e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       force       scalar {1e4}           Load force
%       l           scalar {0.2}           Length of bracket
%       t           scalar {0.02}          Thickness of bracket
%       r           scalar {0.08}          Radius of bracket hole
%       hmax        scalar {0.01}          Max grid cell size
%       sfun        string {sflag2}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        200e9;
  'nu',       0.3;
  'force',    1e4;
  'l',        0.2;
  't',        0.02;
  'r',        0.08-1e-5;
  'hmax',     0.01;
  'sfun',     'sflag2';
  'iplot',    1;
  'tol',      0.05;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.sdim = { 'x' 'y' 'z' };   % Coordinate names.
gobj1    = gobj_block( 0, opt.t, 0, opt.l, 0, opt.l, 'B1' );
gobj2    = gobj_block( 0, opt.l, 0, opt.l, (opt.l-opt.t)/2, (opt.l+opt.t)/2, 'B2' );
gobj3    = gobj_cylinder( [opt.l/2 opt.l/2 opt.l/2-opt.t], opt.r, 2*opt.t, 3, 'C1' );
fea.geom.objects = { gobj1 gobj2 gobj3 };
fea              = geom_apply_formula( fea, 'B1+B2-C1' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid, 'intb', false );


% Problem definition.
fea = addphys(fea,@linearelasticity);
fea.phys.el.eqn.coef{1,end} = { opt.nu };
fea.phys.el.eqn.coef{2,end} = { opt.E  };
fea.phys.el.sfun            = { opt.sfun opt.sfun opt.sfun };


% Boundary conditions.
dtol = 2e-2;
fixbdr   = findbdr( fea, ['x<=',num2str(dtol*1e-1)] );    % Right boundary number.
forcebdr = findbdr( fea, ['x>=',num2str(opt.l-dtol)] );   % Left boundary number.

% Fix right boundary (set zero Dirichlet BCs).
n_bdr  = max(fea.grid.b(3,:));        % Number of boundaries.
bctype = num2cell( zeros(3,n_bdr) );  % First set homogenous Neumann BCs everywhere.
[bctype{:,fixbdr}] = deal( 1 );       % Set Dirchlet BCs for right boundary.
fea.phys.el.bdr.coef{1,5} = bctype;

% Apply negative z-load to left outer boundary.
bccoef = num2cell( zeros(3,n_bdr) );
bccoef{3,forcebdr} = -opt.force;
fea.phys.el.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
warning('off')
fea.sol.u = solvestat( fea, 'fid', fid );
warning('on')


% Postprocessing.
if ( opt.iplot>0 )
  DSCALE = 5000;

  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'u' )
  title( 'x-displacement' )
  view([30 20])

  subplot(2,2,3)
  postplot( fea, 'surfexpr', 'v' )
  title( 'y-displacement' )
  view([30 20])

  subplot(2,2,2)
  postplot( fea, 'surfexpr', 'w' )
  title( 'z-displacement' )
  view([30 20])

  subplot(2,2,4)
  dp = zeros(size(fea.grid.p));
  for i=1:3
    dp(i,:) = DSCALE*evalexpr( fea.dvar{i}, fea.grid.p, fea );
  end
  fea_disp.grid   = fea.grid;
  fea_disp.grid.p = fea_disp.grid.p + dp;
  plotgrid( fea, 'facecolor', [.95 .95 .95], 'edgecolor', [.8 .8 1], ...
                 'selcells', selcells(fea,['x>',num2str(1.5*opt.t)]) )
  hold on
  plotgrid( fea_disp )
  title(['Displacement plot (at ',num2str(DSCALE),' times scale)'])
  view([30 20])

end


% Error check.
xdisp = fea.sol.u(fea.eqn.dofm{1}(:));
ydisp = fea.sol.u(fea.eqn.dofm{2}(:)+fea.eqn.ndof(1));
zdisp = fea.sol.u(fea.eqn.dofm{3}(:)+sum(fea.eqn.ndof(1:2)));
xdisp = [ min(xdisp) max(xdisp) ];
ydisp = [ min(ydisp) max(ydisp) ];
zdisp = [ min(zdisp) max(zdisp) ];
xdisp_ref = [ -8.971e-7 8.971e-7 ];
ydisp_ref = [ -9.556e-8 9.605e-8 ];
zdisp_ref = [ -1.09e-5  1.065e-8 ];
svm_ref = [ 6.252756 1.779065e6 ];
out.xdisp = xdisp;
out.ydisp = ydisp;
out.zdisp = zdisp;
out.err   = abs(zdisp(1) - zdisp_ref(1))/abs(zdisp_ref(1));
out.pass  = out.err<opt.tol;


if ( nargout==0 )
  clear fea out
end
