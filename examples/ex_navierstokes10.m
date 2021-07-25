function [ fea, out ] = ex_navierstokes10( varargin )
%EX_NAVIERSTOKES10 3D Example for stationary flow in a pipe.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES10( VARARGIN ) Sets up and solves
%   stationary and laminar 3D flow in a circular pipe. The inflow
%   profile is constant and the outflow should assume an offset
%   parabolic profile. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.01}          Molecular/dynamic viscosity
%       uin         scalar {0.3}           Magnitude of inlet velocity
%       R           scalar {0.5}           Channel radius
%       sf_u        string {sf_hex_Q1nc}   Shape function for velocity
%       sf_p        string {sf_disc0}      Shape function for pressure
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2 or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = {   ...
  'rho',      1;
  'miu',      1e-2;
  'uin',      0.3;
  'R',        0.5;
  'sf_u',     'sf_hex_Q1nc';
  'sf_p',     'sf_disc0';
  'tol',      0.25;
  'solver',   '';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x', 'y', 'z' };
fea.geom.objects = { gobj_cylinder([0 0 0],opt.R,3,1) };
fea.grid = cylgrid(4,4,20,opt.R,3,[0;0;0],1);


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_u opt.sf_p };
if( any(strcmp(opt.solver,{'openfoam','su2'})) )
  [fea.phys.ns.sfun{:}] = deal('sflag1');
end


% Boundary conditions.
fea.phys.ns.bdr.sel(5) = 2;
fea.phys.ns.bdr.sel(6) = 4;
fea.phys.ns.bdr.coef{2,end}{1,5} = opt.uin;
fea.phys.ns.prop.artstab.iupw = 4;


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid, 'nproc', 1 );
  fid = logfid;
else
  fea.sol.u = solvestat( fea, 'fid', fid );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'sliceexpr', 'sqrt(u^2+v^2+w^2)' )
end


% Error checking.
n = 15;
y = linspace(0.05,0.95,n)' - 0.5;
p = repmat([3 0 0]',1,n);
p(2,:) = y;
u = evalexpr( 'u', p, fea );
u_ref = 2*opt.uin*(1-(y/opt.R).^2);
out.err  = mean(abs(u-u_ref)./u_ref);
out.pass = out.err < opt.tol;


if ( nargout==0 )
  clear fea out
end
