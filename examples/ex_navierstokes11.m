function [ fea, out ] = ex_navierstokes11( varargin )
%EX_NAVIERSTOKES11 3D Example flow in a cubcic cavity.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES11( VARARGIN ) Sets up and solves stationary
%   and laminar 3D flow in a cubic cavity. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {1}             Molecular/dynamic viscosity
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
  'miu',      1;
  'uin',      1;
  'sf_u',     'sf_hex_Q1nc';
  'sf_p',     'sf_disc0';
  'tol',      0.2;
  'solver',   '';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x', 'y', 'z' };
fea.geom.objects = { gobj_block() };
fea.grid = blockgrid(8);
if( strcmp(opt.solver,'openfoam') )
  fea.grid = gridgen( fea, 'hmax', 0.175, 'fid', fid );
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
fea.phys.ns.bdr.sel(6) = 2;
fea.phys.ns.bdr.coef{2,end}{1,6} = opt.uin;


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
u = evalexpr( 'u', [0.5;0.5;0.5], fea );
out.err  = [u-(-0.21)]/0.21;
out.pass = out.err < opt.tol;


if ( nargout==0 )
  clear fea out
end
