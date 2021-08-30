function [ fea, out ] = ex_navierstokes7( varargin )
%EX_NAVIERSTOKES7 3D Example for incompressible stationary flow in a curved pipe.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES7( VARARGIN ) Sets up and solves stationary
%   flow in a curved circular channel. The inflow profile is constant and the outflow
%   should assume an offset parabolic profile. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.001}         Molecular/dynamic viscosity
%       umax        scalar {0.3}           Maximum magnitude of inlet velocity
%       h           scalar {0.5}           Channel radius
%       l           scalar {2.5}           Channel length
%       ilev        scalar {1}             Grid refinement level
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2 or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'rho',      1;
  'miu',      1e-3;
  'umax',     0.3;
  'h',        0.5;
  'l',        0.5;
  'ilev',     1;
  'sf_u',     'sflag2';
  'sf_p',     'sflag1';
  'solver',   '';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Grid generation.
fea.sdim = { 'x' 'y' 'z' };   % Coordinate names.
fea.grid = gridrevolve( circgrid( 4*opt.ilev, 3*opt.ilev, opt.h/2 ), 15*opt.ilev, opt.l, 1/4 );


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
fea.phys.ns.bdr.coef{2,end}{2,5} = -opt.umax;


% Parse and solve problem.
fea       = parsephys(fea);
fea       = parseprob(fea);
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  jac.form  = { [1;1] [1;1] [1;1] []; [1;1] [1;1] [1;1] [];  [1;1] [1;1] [1;1] []; [] [] [] [] };
  jac.coef  = { 'rho_ns*ux' 'rho_ns*uy' 'rho_ns*uz' []; 'rho_ns*vx' 'rho_ns*vy' 'rho_ns*vz' []; 'rho_ns*wx' 'rho_ns*wy' 'rho_ns*wz' []; [] [] [] [] };
  fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)' )
  view( 130, 30 )
end


% Error checking.
out.flow_in  = pi*(opt.h/2)^2*opt.umax;
out.flow_out = intbdr( 'sqrt(u^2+v^2+w^2)', fea, 5 );
out.rerr = abs(out.flow_out-out.flow_in)/out.flow_in;
out.pass = out.rerr < 0.15;


if ( nargout==0 )
  clear fea out
end
