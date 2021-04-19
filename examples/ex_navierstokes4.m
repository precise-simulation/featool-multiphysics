function [ fea, out ] = ex_navierstokes4( varargin )
%EX_NAVIERSTOKES4 2D Example for incompressible flow over a backwards facing step.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES4( VARARGIN ) Stationary flow over a backwards
%   facing step. References:
%
%   [1] P.M. Gresho and R.L. Sani, Incompressible Flow and the Finite Element Method,
%       Volume 1 & 2, John Wiley & Sons, New York, 2000.
%
%   [2] A. Rose and B. Simpson: “Laminar, Constant-Temperature Flow Over a Backward
%       Facing Step,” 1st NAFEMS Workbook of CFD Examples, Glasgow, UK, 2000.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       re          scalar {389}           Reynolds number
%       h           scalar {1}             Channel height
%       y           scalar {0.485}         Step height (fraction of channel height)
%       lc          scalar {7.92}          Channel length (fraction of channel height)
%       li          scalar {1.98}          Inlet length (fraction of channel height)
%       hmax        scalar {0.1}           Max grid cell size
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2 or default solver
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  're',       389;
  'h',        1;
  'y',        0.0049/0.0101;
  'lc',       0.08/0.0101;
  'li',       0.02/0.0101;
  'hmax',     0.1;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'iphys',    1;
  'solver',   '';
  'iplot',    1;
  'tol',      0.2;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid parameters.
h         = opt.h;       % Height of channel.
y         = opt.y;       % Height of step.
lc        = opt.lc;      % Length of channel.
li        = opt.li;      % Length of inlet.
% Model parameters.
umax      = 1;    % Maximum magnitude of inlet velocity.
rho       = 1;           % Density.
miu       = umax*2/3*h/opt.re;    % Molecular/dynamic viscosity.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Geometry definition.
vert = [ -li*h      lc*h lc*h    0 0 -li*h;   ...
         (1-y)*h (1-y)*h -y*h -y*h 0     0];
gobj = gobj_polygon( vert' );
fea.geom.objects = { gobj };
fea.sdim = { 'x' 'y' };    % Coordinate names.


% Grid generation.
fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
n_bdr    = max(fea.grid.b(3,:));        % Number of boundaries.


% Boundary conditions.
dtol      = opt.hmax;
i_inflow  = findbdr( fea, ['x<',num2str(-li*h+dtol)] );   % Inflow boundary number.
i_outflow = findbdr( fea, ['x>',num2str( lc*h-dtol)] );   % Outflow boundary number.
s_inflow  = ['4*',num2str(umax),'*(y*(',num2str((1-y)*h),'-y))/',num2str((1-y)*h),'^2'];   % Definition of inflow profile.
u_init    = ['4*',num2str(umax),'*(y*(',num2str((1-y)*h),'-y))/',num2str((1-y)*h),'^2*(y>0)'];


% Problem definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
  fea.phys.ns.eqn.coef{1,end} = { rho };
  fea.phys.ns.eqn.coef{2,end} = { miu };
  if( ~strcmp(opt.solver,'openfoam') )
    fea.phys.ns.eqn.coef{5,end} = { u_init };
  end
  fea.phys.ns.sfun = { sf_u sf_u sf_p };     % Set shape functions.

  fea.phys.ns.bdr.sel(i_inflow)  = 2;
  fea.phys.ns.bdr.sel(i_outflow) = 4;
  fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;   % Set inflow profile.
  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u'  'v'  'p'  };       % Dependent variable name.
  fea.sfun  = { sf_u sf_u sf_p };       % Shape function.

  % Define equation system.
  cvelx = [num2str(rho),'*',fea.dvar{1}];   % Convection velocity in x-direction.
  cvely = [num2str(rho),'*',fea.dvar{2}];   % Convection velocity in y-direction.
  fea.eqn.a.form = { [2 3 2 3;2 3 1 1]       [2;3]                   [1;2];
                     [3;2]                   [2 3 2 3;2 3 1 1]       [1;3];
                     [2;1]                   [3;1]                   []   };
  fea.eqn.a.coef = { {2*miu miu cvelx cvely}  miu                    -1;
                      miu                    {miu 2*miu cvelx cvely} -1;
                      1                       1                      [] };
  fea.eqn.f.form = { 1 1 1 };
  fea.eqn.f.coef = { 0 0 0 };


  % Define boundary conditions.
  fea.bdr.d = cell(3,n_bdr);
 [fea.bdr.d{1:2,:}]         = deal( 0);

  fea.bdr.d{1,i_inflow}     = s_inflow;

 [fea.bdr.d{:,i_outflow  }] = deal([]);
  % fea.bdr.d{end,i_outflow}  = 0;   % Set pressure to zero on outflow boundary.

  fea.bdr.n = cell(3,n_bdr);
end


% Parse and solve problem.
fea = parseprob(fea);             % Check and parse problem struct.
if( opt.iphys==1 && strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'init', [s_inflow,0,0], 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'fenics') )
  fea = fenics( fea );
else
  fea.sol.u = solvestat(fea,'fid',fid,'maxnit',50,'nlrlx',1,'tolchg',1e-3);   % Call to stationary solver.
end


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr',s_velm,'evaltype','exact','isoexpr',s_velm)
  title('Velocity field')
  subplot(3,1,2)
  postplot(fea,'surfexpr','p','evaltype','exact')
  title('Pressure')
  subplot(3,1,3)
  h = postplot(fea,'surfexpr',['(u<-eps)*x/',num2str(y)]);
  title('Separation length')
end


% Error checking.
s_expr = ['(u<-eps)*x/',num2str(y)];
[~,slen] = minmaxsubd( s_expr, fea );
if( ~isempty(fid) )
  fprintf(fid,'\nRecirculation zone length: %3f (Ref: 7.93)\n\n',slen)
  fprintf(fid,'\n\n')
end

out.slen = [slen 7.93];
out.err  = abs(diff(out.slen))/out.slen(end);
out.pass = out.err<opt.tol;
if ( nargout==0 )
  clear fea out
end
