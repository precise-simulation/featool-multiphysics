function [ fea, out ] = ex_navierstokes12( varargin )
%EX_NAVIERSTOKES12 3D Example flow over a backwards facing step
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES12( VARARGIN ) Sets up and solves stationary and
%   laminar 3D flow over a backwards facing step. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {2/3/389}       Molecular/dynamic viscosity
%       uin         scalar {1}             Magnitude of inlet velocity
%       sf_u        string {sf_hex_Q1nc}   Shape function for velocity
%       sf_p        string {sf_disc0}      Shape function for pressure
%       solver      string 'openfoam'/{''} Use OpenFOAM or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = {   ...
  'rho',      1;
  'miu',      2/3/389;
  'uin',      1;
  'igrid',    1;
  'sf_u',     'sf_hex_Q1nc';
  'sf_p',     'sf_disc0';
  'tol',      0.1;
  'solver',   '';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry.
fea.sdim = { 'x', 'y', 'z' };
gobj1 = gobj_block( -1.9802, 7.9208, 0, 1, 0, 1, 'B1' );
gobj2 = gobj_block( -1.9802, 0, 0, 1, 0, 0.4851, 'B2' );
fea.geom.objects = { gobj1 gobj2 };
fea = geom_apply_formula( fea, 'B1-B2' );


% Grid generation.
if( opt.igrid==1 )
  n = 4;
  fea.grid = rectgrid(10*n,n,[-1.9802, 7.9208;0 1]);
  fea.grid = delcells( fea.grid, selcells(fea.grid,'(y<=0.5).*(x<=0)') );
  ix = find( abs(fea.grid.p(2,:)-0.5)<=sqrt(eps) );
  fea.grid.p(2,ix) = 0.4851;
  fea.grid = gridextrude( fea.grid, n, 1, -2 );
  fea.grid.p(2,:) = fea.grid.p(2,:) + 1;
  fea.grid = assign_bdr( fea.grid, fea.geom, [] );
  [~,ix] = findbdr( fea , 'z>=1-sqrt(eps)', 0 );
  fea.grid.b(3,ix) = 6;
  fea.grid = gridrefine( fea.grid, fid );
else
  fea.grid = gridgen( fea, 'hmax', 0.4, 'fid', fid );
  fea.grid = gridsmooth( tet2hex( fea.grid ), 5 );
end


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_u opt.sf_p };
fea.phys.ns.prop.artstab.iupw = 4;
if( any(strcmp(opt.solver,{'openfoam','su2'})) )
  [fea.phys.ns.sfun{:}] = deal('sflag1');
end


% Boundary conditions.
fea.phys.ns.bdr.sel(5) = 2;
fea.phys.ns.bdr.sel(3) = 4;
fea.phys.ns.bdr.coef{2,end}{1,5} = ['4*',num2str(opt.uin),'*(z-0.4851)*(1-z)/(1-0.4851)^2'];


% Parse and solve problem.
fea  = parsephys( fea );
fea  = parseprob( fea );
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  fea.sol.u = solvestat( fea, 'fid', fid );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'sliceexpr', 'sqrt(u^2+v^2+w^2)' )
end


% Error checking.
[~,slen] = minmaxsubd( '(u<-eps)*x/0.4851*(z<0.4851)*(z>0)*(y<0.51)*(y>0.49)', fea );
out.err  = abs( slen - 3 )/3;
out.pass = out.err < opt.tol;


if ( nargout==0 )
  clear fea out
end
