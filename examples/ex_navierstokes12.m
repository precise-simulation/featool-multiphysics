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
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       solver      string 'openfoam'/{''} Use OpenFOAM or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES4

% Copyright 2013-2023 Precise Simulation, Ltd.


cOptDef = {   ...
  'rho',      1;
  'miu',      2/3/389;
  'uin',      1;
  'igrid',    1;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'solver',   '';
  'iplot',    1;
  'tol',      0.55;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry.
h_step = 0.0049/0.0101;
l_inlet = 0.02/0.0101;
l_channel = 0.08/0.0101;
fea.sdim = { 'x', 'y', 'z' };
gobj1 = gobj_block( -l_inlet, l_channel, -0.5, 0.5, -h_step, 1-h_step, 'B1' );
gobj2 = gobj_block( -l_inlet, 0, -0.5, 0.5, -h_step, 0, 'B2' );
fea.geom.objects = { gobj1, gobj2 };
fea = geom_apply_formula( fea, 'B1-B2' );


% Grid generation.
if( opt.igrid>=1 )
  n = 4;
  fea.grid = rectgrid(10*n,n,[-l_inlet, l_channel; -0.5, 0.5]);
  fea.grid = delcells( fea.grid, selcells(fea.grid,'(y<=0).*(x<=0)') );
  ix = find( fea.grid.p(2,:) <= -0.5 + sqrt(eps) );
  fea.grid.p(2,ix) = -h_step;
  ix = find( fea.grid.p(2,:) >= 0.5 - sqrt(eps) );
  fea.grid.p(2,ix) = 1-h_step;

  fea.grid = gridextrude( fea.grid, n, 1, -2 );
  fea.grid.p(2,:) = fea.grid.p(2,:) + 0.5;
  fea.grid = assign_bdr( fea.grid, fea.geom );
  for i=1:opt.igrid
    fea.grid = gridrefine( fea.grid, fid );
  end
else
  fea.grid = gridgen( fea, 'hmax', 0.1, 'fid', fid );
  % fea.grid = gridsmooth( tet2hex( fea.grid ), 5 );
end


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_u opt.sf_p };
% fea.phys.ns.prop.artstab.iupw = 4;
if( any(strcmp(opt.solver,{'openfoam','su2'})) )
  [fea.phys.ns.sfun{:}] = deal('sflag1');
end


% Boundary conditions.
i_inflow  = findbdr( fea, ['x<',num2str(-l_inlet+1e-3)] );   % Inflow boundary number.
i_outflow = findbdr( fea, ['x>',num2str( l_channel-1e-3)] );   % Outflow boundary number.
% s_inflow  = ['4*',num2str(umax),'*(y*(',num2str((1-y)*h),'-y))/',num2str((1-y)*h),'^2'];   % Definition of inflow profile.
s_inflow  = ['4*',num2str(opt.uin),'*(z*(',num2str(1-h_step),'-z))/(1-',num2str(1-h_step),')^2'];
u_init    = [s_inflow,'*(z>0)'];
fea.phys.ns.bdr.sel(i_inflow) = 2;
fea.phys.ns.bdr.sel(i_outflow) = 4;
fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;
if( ~strcmp(opt.solver,'openfoam') )
  fea.phys.ns.eqn.coef{6,end} = { u_init };
end


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
  fea.sol.u = solvestat( fea, 'maxnit',50, 'nlrlx',1, 'tolchg',1e-3, 'fid', fid );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'sliceexpr', 'sqrt(u^2+v^2+w^2)' )
end


% Error checking.
[~,slen] = minmaxsubd( ['(u<-eps)*x/',num2str(h_step),'*(z<0)*(y<0.01)*(y>-0.01)'], fea );
if( ~isempty(fid) )
  fprintf(fid,'\nRecirculation zone length: %3f (Ref: 7.93)\n\n',slen)
  fprintf(fid,'\n\n')
end

out.slen = [slen, 7.93];
out.err  = abs(diff(out.slen))/out.slen(end);
out.pass = out.err<opt.tol;

if ( nargout==0 )
  clear fea out
end
