function [ fea, out ] = ex_navierstokes8( varargin )
%EX_NAVIERSTOKES8 2D Example for axisymmetric incompressible stationary flow in a constricted circular pipe.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES8( VARARGIN ) Sets up and solves stationary
%   axisymmetric Poiseuille flow in a constricted circular pipe. The inflow
%   profile is constant and the outflow should assume a parabolic profile.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {1}             Molecular/dynamic viscosity
%       uin         scalar {1}             Inflow velocity (constant/mean)
%       r           scalar {1}             Channel radius
%       l           scalar {3}             Channel length
%       hmax        scalar {0.1}           Max grid cell size
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
%   See also EX_NAVIERSTOKES8B

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'rho',      1;
  'miu',      1;
  'uin',      1;
  'r',        1;
  'l',        3;
  'hmax',     0.05;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'solver',   '';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
r = opt.r;   % Pipe radius.
l = opt.l;   % Pipe length.
gobj1 = gobj_rectangle( 0, r,   0,     l*2/3, 'R1' );
gobj2 = gobj_rectangle( 0, r/2, l*2/3, l,     'R2' );
gobj3 = gobj_circle( [r l*2/3], r/2, 'C1' );
fea.geom.objects = { gobj1 gobj2 gobj3 };
fea = geom_apply_formula( fea, 'R1+R2-C1' );
fea.sdim = { 'r' 'z' };   % Space coordinate names.


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );


% Boundary specifications.
i_inflow   = 1;       % Inflow boundary number.
i_outflow  = 5;       % Outflow boundary number.
i_symmetry = [3 6];   % Symmetry boundary numbers.


% Problem definition.
fea = addphys( fea, {@navierstokes 1} );   % Add Navier-Stokes equations physics mode.
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };


% Boundary conditions.
dtol = 1e-3;
i_in  = findbdr( fea, ['z<=',num2str(dtol)] );
i_out = findbdr( fea, ['z>=',num2str(3-dtol)] );
i_sym = findbdr( fea, ['r<=',num2str(dtol)] );
fea.phys.ns.bdr.sel(i_in)   = 2;
fea.phys.ns.bdr.sel(i_out)  = 3;
fea.phys.ns.bdr.coef{2,end}{2,i_in} = opt.uin;
fea.phys.ns.bdr.sel(i_sym) = 5;
fea = parsephys(fea);  % Parse physics mode.


% Parse and solve problem.
fea = parseprob(fea);   % Check and parse problem struct.
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  jac.form  = {[1;1] [1;1] [];[1;1] [1;1] []; [] [] []};
  jac.coef  = {'r*rho_ns*ur' 'r*rho_ns*uz' []; 'r*rho_ns*wr' 'r*rho_ns*wz' []; [] [] []};
  fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac );
end


% Error checking.
r = linspace( 0, opt.r/2, 20 );
z = 0.9*opt.l*ones( 1, 20 );
U = evalexpr( 'sqrt(u^2+w^2)', [r;z], fea )';
u_fac = 4;   % Due to contraction to 1/2 radius.
U_ref = 2*opt.uin*u_fac*( 1 - ( r/(opt.r/2) ).^2 );
err = sqrt( sum((U-U_ref).^2)/sum(U_ref.^2) );


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,3,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+w^2)', 'arrowexpr', {'u' 'w'} )
  hold on
  plot( r, z, 'k--' )
  title( 'Velocity field' )

  subplot(1,3,2)
  postplot( fea, 'surfexpr', 'p', 'evaltype', 'exact', 'isoexpr', 'p' )
  title( 'Pressure' )

  subplot(1,3,3)
  plot( r, U,     'b-' )
  hold on
  plot( r, U_ref, 'r-' )
  title('Velcity profile at z=0.9*l')
  xlabel( 'Radius' )
  legend( 'Computed', 'Reference', 'location', 'south' )
end


out.err  = err;
out.pass = err<0.05;
if( nargout==0 )
  clear fea out
end
