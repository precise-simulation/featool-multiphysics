function [ fea, out ] = ex_navierstokes17( varargin )
%EX_NAVIERSTOKES17 2D Example for incompressible turbulent flow in a channel.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES17( VARARGIN ) Sets up and solves
%   stationary turbulent flow at Re = 42800 in a between two flat
%   parallel plates. The inflow profile is constant and the outflow
%   should assume a fully developed turbulent profile (Ref. Laufer,
%   J. 1950). Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar {3}             Grid refinement level
%                                          (>0=quadrilaterals, <0=triangles)
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.

cOptDef = { ...
            'igrid',    3;
            'solver',   '';
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model parameters.
Re        = 42800;
rho       = 1;           % Density.
miu       = 1/Re;        % Molecular/dynamic viscosity.
uin       = 1;
% Geometry and grid parameters.
h         = 1;           % Height of rectangular domain.
l         = 5;           % Length of rectangular domain.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Geometry definition.
gobj = gobj_rectangle( 0, l, -h/2, h/2 );
fea.geom.objects = { gobj };
fea.sdim = { 'x', 'y' };   % Coordinate names.


% Grid generation.
fea.grid = rectgrid( linspace(0,5,20), [-0.5 -0.45 -0.35 -0.2 0 0.2 0.35 0.45 0.5] );
for i=1:abs(opt.igrid)-1
  fea.grid = gridrefine( fea.grid, fid );
end
if( opt.igrid<0 )
  fea.grid = quad2tri( fea.grid );
end


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
fea.phys.ns.eqn.coef{5,end} = { uin };
if( strcmp(opt.solver,'openfoam') )
  fea.phys.ns.sfun = { 'sflag1', 'sflag1', 'sflag1' };
else
  fea.phys.ns.sfun = { sf_u sf_u sf_p };           % Set shape functions.
end

% Boundary conditions.
i_inflow  = 4;
i_outflow = 2;
fea.phys.ns.bdr.sel(i_inflow)  = 2;
fea.phys.ns.bdr.sel(i_outflow) = 3;
fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = uin;     % Set inflow velocity.

fea.phys.ns.prop.turb.model = 'algebraic';
if( strcmp(opt.solver,'openfoam') )
  fea.phys.ns.prop.turb.model = 'SpalartAllmaras';
end
fea = parsephys(fea);   % Check and parse physics modes.


% Parse and solve problem.
fea = parseprob(fea);             % Check and parse problem struct.
if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  mxnit = 25;
  nlrlx = 0.5;
  told  = 1e-3;
  tolc  = 1e-4;
  fea.sol.u = solvestat( fea, 'fid', fid, 'maxnit', mxnit, 'nlrlx', nlrlx, 'tolchg', tolc, 'toldef', told );
end


% Postprocessing.
y_ref = [0,0.09830508474576277,0.2101694915254238,0.31186440677966104,0.3661016949152544,0.4508474576271188,0.5491525423728814,0.6271186440677967,0.7016949152542376,0.7864406779661017,0.8440677966101697,0.898305084745763,0.9593220338983053,0.9796610169491529,0.9864406779661019,1]/2;
u_ref = [0.9972834919897843,0.9886928256326912,0.9780357557464595,0.9674135128859997,0.9487346180636177,0.9278964476433715,0.9029022521476668,0.8718133271418624,0.8469003947062924,0.7972951009983751,0.7457278848386352,0.6962270722080338,0.5357441374506622,0.32608544230322756,0.19455537497097786,0]';
n = length(y_ref);
x = 0.96*l*ones(1,n);
u = evalexpr( 'u', [x;y_ref], fea );
if( opt.iplot>0 )
  figure
  subplot(2,1,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'isoexpr', 'sqrt(u^2+v^2)' )
  title('Velocity field')
  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'p' )
  title('Pressure')

  figure
  plot( u/max(u), y_ref, 'b.-' )
  hold on
  plot( u_ref, y_ref, 'ro--' )
  ylabel('y')
  xlabel('u/u_{max}')
  grid on
  legend('Computed','Laufer (1950)','location','southwest')
end


% Error checking.
err = sqrt(sum((u-u_ref).^2)/sum(u_ref.^2));

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.pass = err<0.26;
if ( nargout==0 )
  clear fea out
end
