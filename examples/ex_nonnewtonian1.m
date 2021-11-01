function [ fea, out ] = ex_nonnewtonian1( varargin )
%EX_NONNEWTONIAN1 2D Example for non-Newtonian flow in a channel.
%
%   [ FEA, OUT ] = EX_NONNEWTONIAN1( VARARGIN ) Sets up and solves
%   stationary non-Newtonian flow in a rectangular channel. The fluid
%   is governed by a power-law model such that the viscosity can be
%   expressed as mu_eff = mu0*K*D^(n-1). A constant pressure drop
%   is prescribed resulting in an analytical solution for the velocity
%   u(y) = n/(n+1)*(1/(mu0*K)*dp)^(1/n)*((h/2)^((n+1)/n)-abs(h/2-y)^((n+1)/n)).
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       mu0         scalar {1.29684}       Newtonian viscosity
%       dp          scalar {1}             Pressure drop
%       n           scalar {0.25}          Power-law exponent
%       K           scalar {1}             Power-law constant
%       h           scalar {1}             Channel height
%       l           scalar {0.2}           Channel length
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Max grid cell size
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
  'rho',      1;
  'mu0',      1.29684;
  'K',        1;
  'n' ,       0.25;
  'dp' ,      1;
  'h',        1;
  'l',        2;
  'igrid',    1;
  'hmax',     0.05;
  'sf_u',     'sflag2';
  'sf_p',     'sflag1';
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid parameters.
h         = opt.h;       % Height of rectangular domain.
l         = opt.l;       % Length of rectangular domain.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Geometry definition.
gobj = gobj_rectangle( 0, l, 0, h );
fea.geom.objects = { gobj };
fea.sdim = { 'x' 'y' };   % Coordinate names.


% Grid generation.
if ( opt.igrid==1 )
  fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
else
  fea.grid = rectgrid(round(l/opt.hmax),round(h/opt.hmax),[0 l;0 h]);
  if( opt.igrid<0 )
    fea.grid = quad2tri( fea.grid );
  end
end
n_bdr = max(fea.grid.b(3,:));           % Number of boundaries.


% Boundary conditions.
dtol      = opt.hmax;
i_inflow  = findbdr( fea, ['x<',num2str(dtol)] );     % Inflow boundary number.
i_outflow = findbdr( fea, ['x>',num2str(l-dtol)] );   % Outflow boundary number.


% Problem definition.
fea = addphys(fea,@nonnewtonian);      % Add non-Newtonian flow physics mode.
fea.phys.nn.eqn.coef{1,end} = { opt.rho };
fea.phys.nn.eqn.coef{2,end} = { 1 };   % Select power-law viscosity model.
fea.phys.nn.eqn.coef{4,end} = { opt.mu0 };
fea.phys.nn.eqn.coef{6,end} = { opt.K };
fea.phys.nn.eqn.coef{7,end} = { opt.n };
fea.phys.nn.eqn.coef{13,end}{1} = [num2str(opt.dp),'*(1-x/',num2str(l),')'];   % Initial pressure.
fea.phys.nn.sfun = { sf_u sf_u sf_p };           % Set shape functions.
fea.phys.nn.bdr.sel([i_inflow,i_outflow]) = 4;
fea.phys.nn.bdr.coef{4,end}{3,i_inflow}   = l/opt.dp;   % Set inflow pressure.


% Parse and solve problem.
fea = parsephys(fea);   % Check and parse physics mode.
fea = parseprob(fea);   % Check and parse problem struct.

[fea.bdr.n{1}{[i_inflow,i_outflow]}] = deal( '-p*nx' );   % Offset natural Neumann BC.
[fea.bdr.d{2}{[i_inflow,i_outflow]}] = deal(0);           % Set v=0 at inlet and outlet.

fea.sol.u = solvestat( fea, 'relchg', false, 'fid', fid, 'maxnit', 50 );   % Call to stationary solver.


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
s_refsol  = [num2str(opt.n/(opt.n+1)*(1/(opt.mu0*opt.K)*opt.dp)^(1/opt.n)),'*(',num2str((h/2)^((opt.n+1)/opt.n)),'-abs(',num2str(h/2),'-y)^',num2str((opt.n+1)/opt.n),')'];
p = [ l/2*ones(1,25); linspace(0,h,25) ];
u = evalexpr( s_velm, p, fea );
u_ref = evalexpr( s_refsol, p, fea );
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr',s_velm,'evaltype','exact')
  title('Velocity field')

  subplot(3,1,2)
  postplot(fea,'surfexpr','miu_nn')
  title('Effective viscosity')

  subplot(3,1,3)
  plot( u, p(2,:) )
  hold on
  plot( u_ref, p(2,:), 'k.' )
  legend( 'Computed solution', 'Reference solution','Location','West')
  title( 'Velocity profile at x=l/2' )
  ylabel( 'y' )
end


% Error checking.
err = sqrt(sum((u-u_ref).^2)/sum(u_ref.^2));
if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.pass = err<0.05;
if ( nargout==0 )
  clear fea out
end
