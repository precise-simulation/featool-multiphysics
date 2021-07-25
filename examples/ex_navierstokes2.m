function [ fea, out ] = ex_navierstokes2( varargin )
%EX_NAVIERSTOKES2 2D Example for incompressible flow in a square cavity.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES2( VARARGIN ) Sets up and solves stationary
%   incompressible flow in a square cavity. References:
%
%   [1] Botella O, Peyret R. Benchmark spectral results on the lid-
%       driven cavity flow. Computers and Fluids 27(4):421–433, 1998.
%
%   [2] Erturk E, Corke TC, Gökcöl C. Numerical solutions of 2-D steady
%       incompressible driven cavity flow at high Reynolds numbers. Int-
%       ernational Journal for Numerical Methods in Fluids 37(6):633–655, 2005.
%
%   [3] Nishida H, Satofuka N. Higher-order solutions of square driven cavity
%       flow using a variable-order multi-grid method. International Journal
%       for Numerical Methods in Engineering 34(2):637–653, 1992.
%
%   [4] Schreiber R, Keller HB. Driven cavity flows by efficient numerical
%       techniques. Journal of Computational Physics 49(2):310–333, 1983.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       re          scalar {1000}          Reynolds number
%       igrid       scalar 0/{1}/2         Cell type (0=quadrilaterals, 1=triangles,
%                                          2=triangles converted from quadrilaterals)
%       hmax        scalar {0.02}          Max grid cell size
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iphys       scalar {1}/0           Use physics mode to define problem (=1)
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2 or default solver
%       iplot       scalar {1}/0           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES2B

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  're',       1000;
  'igrid',    1;
  'hmax',     0.02;
  'sf_u',     'sflag1';
  'sf_p',     'sflag1';
  'iphys',    1;
  'solver',   '';
  'iplot',    1;
  'tol',      0.35;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model parameters.
rho       = 1;           % Density.
umax      = 1;           % Maximum magnitude of inlet velocity.
l         = 1;
miu       = umax*l/opt.re;   % Molecular/dynamic viscosity.
% Grid parameters.
hmax      = opt.hmax;    % Max allowable global element size.
hmaxr     = 2*hmax;      % Max allowable element size on rectangle.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Geometry definition.
gobj = gobj_rectangle( 0, l, 0, l );
gobj.boundaries = gobj.boundaries([4,1,2,3]);
fea.geom.objects = { gobj };
fea.sdim = { 'x' 'y' };   % Coordinate names.


% Grid generation.
switch opt.igrid
  case  1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
  otherwise
    fea.grid = rectgrid(round(l/opt.hmax),round(l/opt.hmax),[0 l;0 l]);
    ib = zeros(1,size(fea.grid.b,2));
    ib(fea.grid.b(3,:)==4) = 1;
    ib(fea.grid.b(3,:)==1) = 2;
    ib(fea.grid.b(3,:)==2) = 3;
    ib(fea.grid.b(3,:)==3) = 4;
    fea.grid.b(3,:) = ib;
    if( opt.igrid~=0 )
      fea.grid = quad2tri(fea.grid,1);
    end
end


% Boundary conditions.
n_bdr    = max(fea.grid.b(3,:));                     % Number of boundaries.
dtol     = opt.hmax;
i_inflow = findbdr( fea, ['y>',num2str(l-dtol)] );   % Inflow (top) boundary.


% Problem definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
  fea.phys.ns.eqn.coef{1,end} = { rho };
  fea.phys.ns.eqn.coef{2,end} = { miu };
  fea.phys.ns.sfun            = { sf_u sf_u sf_p };           % Set shape functions.

  fea.phys.ns.bdr.sel(i_inflow) = 2;
  fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = umax;             % Set inflow profile.
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
  fea.bdr.d{1,i_inflow}     = umax;
  fea.bdr.n = cell(3,n_bdr);

  % Add integral pressure constraint.
  fea.constr.type  = 'intsubd';
  fea.constr.dvar  = 'p';
  fea.constr.index = [];
  fea.constr.expr  = 0';
end


% Parse and solve problem.
fea = parseprob(fea);             % Check and parse problem struct.
if( opt.iphys==1 && strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'tol', 3e-7, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', fid, 'nlrlx', 0.5, 'maxnit', 50 );
else
  fea.sol.u = solvestat( fea, 'fid', fid, 'maxnit', 50 );   % Call to stationary solver.
end


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr',s_velm,'evaltype','exact')
  title('Velocity field')
  subplot(3,1,2)
  postplot(fea,'surfexpr','p','evaltype','exact')
  title('Pressure')
  subplot(3,1,3)
  postplot(fea,'surfexpr','vx-uy','evaltype','exact','isoexpr',s_velm,'isolev',30)
  title('Vorticity')
end


% Error checking.
vort     = evalexpr('vx-uy',[0.53;0.564],fea);
out.aerr = -2.068-vort;
out.err  = abs(out.aerr)/2.068;
out.pass = out.err < opt.tol;
if ( nargout==0 )
  clear fea out
end
