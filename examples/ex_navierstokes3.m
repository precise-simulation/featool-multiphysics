function [ fea, out ] = ex_navierstokes3( varargin )
%EX_NAVIERSTOKES3 2D Example for incompressible stationary flow around a cylinder.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES3( VARARGIN ) Stationary flow around a cylinder. References:
%
%   [1] John V, Matthies G. Higher-order finite element discretizations in a
%       benchmark problem for incompressible flows. International Journal for
%       Numerical Methods in Fluids 2001.
%
%   [2] Nabh G. On higher order methods for the stationary incompressible
%       Navier-Stokes equations. PhD Thesis, Universitaet Heidelberg, 1998.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.001}         Molecular/dynamic viscosity
%       umax        scalar {0.3}           Maximum magnitude of inlet velocity
%       igrid       scalar {2}             Grid type: >0 regular (igrid refinements)
%                                                     <0 unstruc. grid (with hmax=|igrid|)
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       solver      string {}              Solver selection default, openfoam, su2, fenics
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES3B

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
            'rho',      1;
            'miu',      0.001;
            'umax',     0.3;
            'igrid',    2;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iphys',    1;
            'solver',   '';
            'iplot',    1;
            'tol',      [0.05 0.35 0.1];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model parameters.
rho       = opt.rho;     % Density.
miu       = opt.miu;     % Molecular/dynamic viscosity.
umax      = opt.umax;    % Maximum magnitude of inlet velocity.
umean     = 2/3*umax;    % Mean inlet velocity.
% Geometry and grid parameters.
h         = 0.41;        % Height of rectangular domain.
l         = 2.2;         % Length of rectangular domain.
xc        = 0.2;         % x-coordinate of cylinder center.
yc        = 0.2;         % y-coordinate of cylinder center.
diam      = 0.1;         % Diameter of cylinder.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Grid generation.
fea.sdim = { 'x' 'y' };
if( opt.igrid>=1 )
  fea.grid = cylbenchgrid( opt.igrid );
  fea.grid.s(:) = 1;
else
  gobj1 = gobj_rectangle( 0, 2.2, 0, 0.41, 'R1' );
  gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
  fea.geom.objects = { gobj1 gobj2 };
  fea = geom_apply_formula( fea, 'R1-C1' );
  fea.grid = gridgen( fea, 'hmax', abs(opt.igrid), 'fid', fid );
end
if( opt.iphys==1 && strcmp(opt.solver,'fenics') &&  size(fea.grid.c,1)==4 )
  warning( 'Converting quadrilateral mesh to triangular to support FEniCS.' )
  fea.grid = quad2tri( fea.grid );
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Boundary conditions.
dtol      = sqrt(eps)*1e3;
i_inflow  = findbdr( fea, ['x<=',num2str(dtol)] );     % Inflow boundary number.
i_outflow = findbdr( fea, ['x>=',num2str(l-dtol)] );   % Outflow boundary number.
s_inflow  = ['4*',num2str(umax),'*(y*(',num2str(h),'-y))/',num2str(h),'^2'];   % Definition of inflow profile.
i_cyl     = findbdr( fea, ['sqrt((x-',num2str(xc),').^2+(y-',num2str(yc),').^2)<=(',num2str(diam/2+dtol),')'] );    % Cylinder boundary number.


% Problem definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
  fea.phys.ns.eqn.coef{1,end} = { rho };
  fea.phys.ns.eqn.coef{2,end} = { miu };
  fea.phys.ns.sfun            = { sf_u sf_u sf_p };     % Set shape functions.

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
                     miu                     {miu 2*miu cvelx cvely} -1;
                     1                       1                      [] };
  fea.eqn.f.form = { 1 1 1 };
  fea.eqn.f.coef = { 0 0 0 };


  % Define boundary conditions.
  fea.bdr.d = cell(3,n_bdr);
  [fea.bdr.d{1:2,:}]         = deal( 0);

  fea.bdr.d{1,i_inflow}      = s_inflow;

  [fea.bdr.d{:,i_outflow  }] = deal([]);

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
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( opt.iphys==1 && strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', fid );
else
  jac.form  = {[1;1] [1;1] [];[1;1] [1;1] []; [] [] []};
  jac.coef  = {[num2str(rho),'*ux'] [num2str(rho),'*uy'] []; [num2str(rho),'*vx'] [num2str(rho),'*vy'] []; [] [] []};
  fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac, 'nlrlx', '(1+(it>2))/2' );   % Call to stationary solver.
end


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot>0 )
  figure
  subplot(2,1,1)
  postplot( fea, 'surfexpr', s_velm )
  title( 'Velocity field' )
  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'p' )
  title( 'Pressure' )
end


% Calculate benchmark quantities (line integration method).
s_tfx = ['nx*p+',num2str(miu),'*(-2*nx*ux-ny*(uy+vx))'];
s_tfy = ['ny*p+',num2str(miu),'*(-nx*(vx+uy)-2*ny*vy)'];
s_cd  = ['2*(',s_tfx,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
s_cl  = ['2*(',s_tfy,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
i_cub = 10;
c_d1  = intbdr(s_cd,fea,i_cyl,i_cub);
c_l1  = intbdr(s_cl,fea,i_cyl,i_cub);
dp    = evalexpr('p',[0.15 0.25;0.2 0.2],fea);


% Calculate benchmark quantities (volume integration method).
bdrm   = fea.bdr.bdrm{1};
ind_b  = [];
ind_bm = [];
for ii=i_cyl
  ind_b  = [ind_b  find(fea.grid.b(3,:)==ii)];
  ind_bm = [ind_bm find(bdrm(3,:)==ii)];
end
ind_c    = fea.grid.b(1,ind_b);
ind_gdof = bdrm(4,ind_bm);

% Create field 'a' with values one on the cylinder and zero everywhere else.
fea.dvar = [ fea.dvar, {'a'}       ];
fea.sfun = [ fea.sfun, fea.sfun(1) ];
fea      = parseprob(fea);
n_dof    = max(fea.eqn.dofm{1}(:));
u_a      = zeros(n_dof,1);
u_a(ind_gdof) = 1;
fea.sol.u= [fea.sol.u;u_a];
fea.eqn  = struct;
fea.bdr  = struct;
fea      = parseprob(fea);

s_tfx    = ['ax*p+',num2str(miu),'*(-2*ax*ux-ay*(uy+vx))-(u*ux+v*uy)*a'];
s_tfy    = ['ay*p+',num2str(miu),'*(-ax*(vx+uy)-2*ay*vy)-(u*vx+v*vy)*a'];
s_cd     = ['2*(',s_tfx,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
s_cl     = ['2*(',s_tfy,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
c_d2     = intsubd(s_cd,fea,[],[],3);
c_l2     = intsubd(s_cl,fea,[],[],3);


if( ~isempty(fid) )
  fprintf(fid,'\n\nBenchmark quantities:\n\n')

  fprintf(fid,'Drag coefficient,    cd = %6f (l), %6f (v) (Ref: 5.579535)\n',c_d1,c_d2)
  fprintf(fid,'Lift coefficient,    cl = %6f (l), %6f (v) (Ref: 0.010619)\n',c_l1,c_l2)
  fprintf(fid,'Pressure,            dp = %6f (Ref: 0.117520)\n',dp(1)-dp(2))
end


% Error checking.
out.cd   = [c_d1 c_d2];
out.cl   = [c_l1 c_l2];
out.dp   = dp(1)-dp(2);
out.err  = [abs(out.cd-5.579535)/5.579535;
            abs(out.cl-0.010619)/0.010619;
            abs(dp(1)-dp(2)-0.117520)/0.117520 0];
out.pass = (out.err(1,2)<opt.tol(1))&&(out.err(2,2)<opt.tol(2))&&(out.err(3,1)<opt.tol(3));
if ( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ grid ] = cylbenchgrid( nlev )
% CYLBENCHGRID Generate 2d quadrilateral grid for the DFG cylinder benchmark.

  ns = 8*2^(nlev-1);
  r  = [0.05 0.06 0.08 0.11 0.15];
  x  = [0.41 0.5 0.7 1 1.4 1.8 2.2];
  for ilev=2:nlev
    r = sort( [ r (r(1:end-1)+r(2:end))/2 ] );
    x = sort( [ x (x(1:end-1)+x(2:end))/2 ] );
  end

  grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
  grid2 = holegrid( ns, 2^(nlev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
  grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
  grid3 = rectgrid( x, ns, [0.41 2.2;0 0.41] );
  grid  = gridmerge( grid3, 4, grid2, 6 );
