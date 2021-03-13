function [ fea, out ] = ex_darcy1( varargin )
%EX_DARCY1 Darcy's law for a porous packed bed reactor.
%
%   [ FEA, OUT ] = EX_DARCY1( VARARGIN ) Example using Darcy's law to study the
%   cross sectional flow in the porous material of a packed bed reactor.
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


cOptDef = { 'hmax',     2e-4;
            'sfun',     'sflag2';
            'iplot',    1;
            'tol',      5e-3;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model constants.
p0  = 1e5;       % Reference pressure.
kap = 1.2e-13;   % Permeability.
eta = 1.6e-5;    % Viscosity.
M   = 16e-3;     % Molar mass.
R   = 8.32;      % Ideal gas constant.
T   = 344;       % Temperature.


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.geom.objects = { gobj_polygon([0     2.3e-3 2.3e-3 2.3e-3 2.3e-3 0    0    0       0;
                                  -6e-3 -6e-3  -5e-3  -4e-3   6e-3   6e-3 5e-3 3.5e-3 -6e-3]') };
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );


% Problem definition.
fea = addphys( fea, @darcyslaw );   % Add Darcy's law physics mode.

% Modify equation kappa -> p*kappa (introduce nonlinearity in equation, not coefficient).
fea.phys.dl.eqn.seqn = strrep(fea.phys.dl.eqn.seqn,'kap_dl','p*kap_dl');
fea.phys.dl.eqn.coef{2,end} = { kap*M/R/T };
fea.phys.dl.eqn.coef{3,end} = { eta };
fea.phys.dl.sfun            = { opt.sfun };

% Set pressure 1.7 times higher on the inflow.
i_inflow  = 7;
i_outflow = 3;
fea.phys.dl.bdr.sel(i_inflow)  = 1;
fea.phys.dl.bdr.sel(i_outflow) = 1;
fea.phys.dl.bdr.coef{1,end}{1,i_inflow}  = p0*1.7;
fea.phys.dl.bdr.coef{1,end}{1,i_outflow} = p0;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );             % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'init', {p0}, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
s_velmag = ['sqrt(((',num2str(-kap/eta),')^2)*(px^2+py^2))'];
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', ['min(0.1,',s_velmag,')'], ...
                 'arrowexpr', {[num2str(-kap/eta),'*px'] [num2str(-kap/eta),'*py']} )
  title( 'Velocity field' )
end


% Error checking.
u_vel = intsubd( s_velmag, fea );
p_in  = intbdr( 'p', fea, i_inflow );
p_out = intbdr( 'p', fea, i_outflow );
out.err = [ abs(u_vel-1.205e-6)/1.205e-6;
            abs(p_in-255)/255 ;
            abs(p_out-100)/100 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
