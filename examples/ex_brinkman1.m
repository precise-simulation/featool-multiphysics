function [ fea, out ] = ex_brinkman1( varargin )
%EX_BRINKMAN1 Coupled Navier-Stokes and Brinkman equations model.
%
%   [ FEA, OUT ] = EX_BRINKMAN1( VARARGIN ) Example using the Navier-Stokes equations coupled
%   with the Brinkman equations for porous media flow.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {2e-4}          Grid size
%       sfun        string {sflag2}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     4e-4;
            'iplot',    1;
            'tol',      0.25;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model constants.
rho = 1.23;     % Density.
miu = 2.1e-4;   % Viscosity.
kap = 3.2e-7;   % Permeability.


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.geom.objects = { gobj_rectangle(-2e-3,0,-8e-3,8e-3,'R1') ...
                     gobj_polygon([0 -6e-3;2e-3 -5e-3;2e-3 4e-3;0 6e-3],'P1') ...
                                    };
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );
fea.grid = gridbdrx( fea.grid );   %  Reconstruct internal boundaries.


% Problem definition.

% Navier-Stokes equations.
fea = addphys( fea, @navierstokes );   % Add Navier-Stokes equations physics mode.

fea.phys.ns.eqn.coef{1,end} = rho;
fea.phys.ns.eqn.coef{2,end} = miu;

fea.phys.ns.bdr.sel(1) = 2;
fea.phys.ns.bdr.coef{2,end}{2,1} = 2e-2;   % Set inflow velocity.
fea.phys.ns.bdr.sel(4) = 4;
fea.phys.ns.bdr.sel(9) = -2;
fea.phys.ns.bdr.coefi{2,end}{1,9} = 'u2';   % Brinkman velocity from shared boundary.
fea.phys.ns.bdr.coefi{2,end}{2,9} = 'v2';

% Brinkman equations.
fea = addphys( fea, @brinkmaneqns );   % Add Brinkman equations physics mode.

fea.phys.br.eqn.coef{1,end} = rho;
fea.phys.br.eqn.coef{2,end} = miu;
fea.phys.br.eqn.coef{3,end} = kap;

fea.phys.br.bdr.sel(9) = -2;
fea.phys.br.bdr.coefi{2,end}{1,9} = 'u';   % NS velocity from shared boundary.
fea.phys.br.bdr.coefi{2,end}{2,9} = 'v';


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'solcomp', {{1 0} {1 0} {1 0} {0 1} {0 1} {0 1} }, 'fid', opt.fid );


% Postprocessing.
x = linspace( -1e-3, 1e-3, 100 );
y = zeros(size(x));
u_ns = evalexpr( 'sqrt(u^2+v^2)', [x;y], fea );
u_br = evalexpr( 'sqrt(u2^2+v2^2)', [x;y], fea );
if( opt.iplot )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'} )
  postplot( fea, 'surfexpr', 'sqrt(u2^2+v2^2)', 'arrowexpr', {'u2' 'v2'} )

  subplot(1,2,2)
  plot( x, u_ns ), hold on
  plot( x, u_br ), hold on
end


% Error checking.
out.err = [ abs(max( u_ns ) - 0.0175)/0.0175, ...
            abs(mean( u_br(~isnan(u_br)) ) - 9e-3)/9e-3];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
