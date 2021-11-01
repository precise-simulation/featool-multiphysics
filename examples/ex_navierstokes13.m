function [ fea, out ] = ex_navierstokes13( varargin )
%EX_NAVIERSTOKES13 3D Example flow past a cylinder
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES13( VARARGIN ) Sets up and solves stationary and
%   laminar 3D flow past a cylinder. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.001}         Molecular/dynamic viscosity
%       uin         scalar {0.45}          Magnitude of inlet velocity
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
  'miu',      0.001;
  'uin',      0.45;
  'nlev',     1;
  'sf_u',     'sf_hex_Q1nc';
  'sf_p',     'sf_disc0';
  'tol',      0.2;
  'solver',   'openfoam';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry.
fea.sdim = { 'x', 'y', 'z' };
gobj1 = gobj_block( 0, 2.5, 0, 0.41, 0, 0.41, 'B1' );
gobj2 = gobj_cylinder( [0.5 0 0.21], 0.05, 0.41, 2, 'C1' );
fea.geom.objects = { gobj1 gobj2 };
fea = geom_apply_formula( fea, 'B1-C1' );


% Grid generation.
ns = 8*2^(opt.nlev-1);
r  = [0.05 0.06 0.08 0.11 0.15];
x  = [0.41 0.5 0.7 1 1.4 1.8 2.2] + 0.3;
for ilev=2:opt.nlev
  r = sort( [ r (r(1:end-1)+r(2:end))/2 ] );
  x = sort( [ x (x(1:end-1)+x(2:end))/2 ] );
end

grid1 = ringgrid( r, 4*ns, [], [], [0.5;0.2] );
grid2 = holegrid( ns, 2^(opt.nlev-1), [0.3 0.71;0 0.41], 0.15, [0.5;0.2] );
grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
grid3 = rectgrid( x, ns, [0.71 2.5;0 0.41] );
fea.grid = gridmerge( grid3, 4, grid2, 6 );
grid4 = rectgrid( 1, ns, [0 0.3;0 0.41] );
fea.grid = gridmerge( fea.grid, 10, grid4, 2 );
fea.grid = gridextrude( fea.grid, 5, 0.41, 2 );
fea.grid.p(3,:) = fea.grid.p(3,:) + 0.41;

% Renumber boundaries to match geometry.
fea.grid.b(3,fea.grid.b(3,:)==6) = -7;
fea.grid.b(3,fea.grid.b(3,:)==5) = -8;
fea.grid.b(3,fea.grid.b(3,:)==4) = -9;
fea.grid.b(3,fea.grid.b(3,:)==7) = -10;
[~,ind] = findbdr( fea, 'x<=sqrt(eps)', false );
fea.grid.b(3,ind) = 5;
[~,ind] = findbdr( fea, 'x>=2.5-sqrt(eps)', false );
fea.grid.b(3,ind) = 3;
[~,ind] = findbdr( fea, 'y<=sqrt(eps)', false );
fea.grid.b(3,ind) = 2;
[~,ind] = findbdr( fea, 'y>=0.41-sqrt(eps)', false );
fea.grid.b(3,ind) = 4;
[~,ind] = findbdr( fea, 'z<=sqrt(eps)', false );
fea.grid.b(3,ind) = 1;
[~,ind] = findbdr( fea, 'z>=0.41-sqrt(eps)', false );
fea.grid.b(3,ind) = 6;
fea.grid.b(3,fea.grid.b(3,:)==-7)  = 7;
fea.grid.b(3,fea.grid.b(3,:)==-8)  = 8;
fea.grid.b(3,fea.grid.b(3,:)==-9)  = 9;
fea.grid.b(3,fea.grid.b(3,:)==-10) = 10;
fea.grid.s(:) = 1;
if( strcmp(opt.solver,'openfoam') )
  fea.grid = gridrefine( fea.grid, fid );
end


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_u opt.sf_p };
if( strcmp(opt.solver,'openfoam') )
  [fea.phys.ns.sfun{:}] = deal('sflag1');
end


% Boundary conditions.
fea.phys.ns.bdr.sel(5) = 2;
fea.phys.ns.bdr.coef{2,end}{1,5} = ['16*',num2str(opt.uin),'*(y*z*(0.41-y)*(0.41-z))/0.41^4'];
fea.phys.ns.bdr.sel(3) = 4;
fea.phys.ns.prop.artstab.iupw = 4;


% Parse and solve problem.
fea  = parsephys( fea );
fea  = parseprob( fea );
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  fea.sol.u = solvestat( fea, 'fid', fid );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'sliceexpr', 'sqrt(u^2+v^2+w^2)' )
end


% Error checking.
s_tfx = ['nx*p+',num2str(opt.miu),'*(-2*nx*ux-nz*(uz+vx))'];
s_cd  = ['-2*(',s_tfx,')/(',num2str(opt.rho),'*',num2str(0.2),'^2*',num2str(0.1*0.41),')'];
c_d1  = abs(intbdr(s_cd,fea,[7:10],10));
dp    = evalexpr('p',[0.45 0.55;0.205 0.205;0.21 0.21],fea);
out.err = [abs(c_d1-6.185333)/6.185333, ...
           abs(dp(1)-dp(2)-0.171007)/0.171007];
out.pass = all(out.err < [0.15 0.6]);


if ( nargout==0 )
  clear fea out
end
