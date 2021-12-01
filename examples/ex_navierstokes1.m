function [ fea, out ] = ex_navierstokes1( varargin )
%EX_NAVIERSTOKES1 2D Example for incompressible stationary flow in a channel.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES1( VARARGIN ) Sets up and solves stationary Poiseuille
%   flow in a rectangular channel. The inflow profile is constant and the outflow
%   should assume a parabolic profile ( u(y)=U_max*4/h^2*y*(h-y) ).
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.001}         Molecular/dynamic viscosity
%       umax        scalar {0.3}           Maximum magnitude of inlet velocity
%       h           scalar {0.5}           Channel height
%       l           scalar {2.5}           Channel length
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.04}          Max grid cell size
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2, FEniCS, or default solver
%       ischeme     scalar {0}             Time stepping scheme (0 = stationary)
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES1B

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'rho',      1;
  'miu',      1e-3;
  'umax',     0.3;
  'h',        0.5;
  'l',        2.5;
  'igrid',    1;
  'hmax',     0.04;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'iphys',    1;
  'solver',   '';
  'ischeme',  0;
  'iplot',    1;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model parameters.
rho       = opt.rho;     % Density.
miu       = opt.miu;     % Molecular/dynamic viscosity.
umax      = opt.umax;    % Maximum magnitude of inlet velocity.
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
s_inflow  = ['2/3*',num2str(umax)];                                            % Definition of inflow profile.
s_refsol  = ['4*',num2str(umax),'*(y*(',num2str(h),'-y))/',num2str(h),'^2'];   % Definition of velocity profile.


% Problem definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
  fea.phys.ns.eqn.coef{1,end} = { rho };
  fea.phys.ns.eqn.coef{2,end} = { miu };
  fea.phys.ns.eqn.coef{5,end} = { s_inflow };
  if( any(strcmp(opt.solver,{'openfoam','su2'})) )
    fea.phys.ns.sfun = { 'sflag1', 'sflag1', 'sflag1' };
  else
    fea.phys.ns.sfun = { sf_u sf_u sf_p };           % Set shape functions.
  end
  fea.phys.ns.bdr.sel(i_inflow)  = 2;
  fea.phys.ns.bdr.sel(i_outflow) = 4;
  fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;         % Set inflow profile.
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
 [fea.bdr.d{1:2,:}]         = deal( 0 );

  fea.bdr.d{1,i_inflow}     = s_inflow;

 [fea.bdr.d{:,i_outflow  }] = deal([]);
  % fea.bdr.d{end,i_outflow}  = 0;   % Set pressure to zero on outflow boundary.

  fea.bdr.n = cell(3,n_bdr);
end


% Parse and solve problem.
fea = parseprob(fea);             % Check and parse problem struct.
if( opt.iphys==1 && strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', fid, 'ischeme', opt.ischeme, 'tmax', 10 );
elseif( opt.iphys==1 && strcmp(opt.solver,'openfoam') )
  if( opt.ischeme==0 )
    dt = 1.0;
    tstop = 1000;
    ddtScheme = 'steadyState';
  elseif( opt.ischeme==1 )
    dt = 0.1;
    tstop = 100;
    ddtScheme = 'backward';
  elseif( opt.ischeme>=2 )
    dt = 0.1;
    tstop = 100;
    ddtScheme = 'CrankNicolson 0.9';
  end
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid, 'ddtScheme', ddtScheme, 'deltaT', dt, 'endTime', tstop );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid, 'ischeme', opt.ischeme, 'tstep', 0.5, 'tmax', 20+30*(opt.ischeme==1) );
  fid = logfid;
else
  if( opt.ischeme==0 )
    jac.form  = {[1;1] [1;1] [];[1;1] [1;1] []; [] [] []};
    jac.coef  = {[num2str(rho),'*ux'] [num2str(rho),'*uy'] []; [num2str(rho),'*vx'] [num2str(rho),'*vy'] []; [] [] []};
    fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac );   % Call to stationary solver.
  else
    fea.sol.u = solvetime( fea, 'fid', fid, 'ischeme', opt.ischeme, 'tmax', 10 );
  end
end
fea.sol.u = fea.sol.u(:,end);


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
s_err  = ['abs(sqrt((',s_refsol,')^2)-(',s_velm,'))'];
s_len  = ['(x>',num2str(3/4*l),')'];
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr',s_velm,'evaltype','exact')
  title('Velocity field')
  subplot(3,1,2)
  postplot(fea,'surfexpr','p','evaltype','exact')
  title('Pressure')
  subplot(3,1,3)
  postplot(fea,'surfexpr',[s_err,'*',s_len],'evaltype','exact')
  title('Error')
end


% Error checking.
if ( size(fea.grid.c,1)==4 )
  xi = [0;0];
else
  xi = [1/3;1/3;1/3];
end
c_ind = find(evalexpr0(s_len,xi,1,1:size(fea.grid.c,2),[],fea))';
err = evalexpr0(s_err,xi,1,c_ind,[],fea);
ref = evalexpr0(['sqrt((',s_refsol,')^2)'],xi,1,c_ind,[],fea);
err = sqrt(sum(err.^2)/sum(ref.^2));


if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.pass = err<0.06;
if ( nargout==0 )
  clear fea out
end
