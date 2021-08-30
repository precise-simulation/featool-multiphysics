function [ data, fea ] = ex_navierstokes3b( varargin )
%EX_NAVIERSTOKES3B 2D Flow around cylinder CFD benchmark
%
%   [ DATA, FEA ] = EX_NAVIERSTOKES3B 2D Flow around cylinder CFD benchmark
%
%   CFD benchmarking script comparing the OpenFOAM and SU2 CFD solvers
%   to the fully coupled and monolithic FEniCS and FEATool
%   Multiphysics solvers, for a 2D stationary and incompressible flow
%   around a cylinder test case.
%
%   Benchmark references:
%
%   [1] John V, Matthies G. Higher-order finite element
%       discretizations in a benchmark problem for incompressible flows.
%       International Journal for Numerical Methods in Fluids 2001.
%
%   [2] Nabh G. On higher order methods for the stationary
%       incompressible Navier-Stokes equations. PhD Thesis,
%       Universitaet Heidelberg, 1998.
%
%   See also EX_NAVIERSTOKES3, RUN_FEATOOL_BENCHMARKS

% Copyright 2013-2021 Precise Simulation, Ltd.


% Define available solvers.
solvers  = { 'FEATool', 'FEniCS', 'OpenFOAM', 'SU2' };
FEATOOL  = 1;
FENICS   = 2;
OPENFOAM = 3;
SU2      = 4;

% Define available grid cases.
grids    = { 'Quad', 'Tri', 'TriU' };
QUAD     = 1;
TRI      = 2;   % Structured triangles.
TRIU     = 3;   % Unstructured triangles.

% Define available FEM shape/basis functions.
D0       = 'sf_disc0';
D1       = 'sf_disc1';
P1       = 'sf_simp_P1';
P2       = 'sf_simp_P2';
Q1       = 'sf_quad_Q1';
Q2       = 'sf_quad_Q2';

% Set up benchmark test cases.
NLEV  = 1:5;   % Select grid levels.
cases = { FEATOOL,   TRI,   {P1,P1},  NLEV;
          FEATOOL,   TRI,   {P2,P1},  NLEV(1:end-1);
          FEATOOL,   QUAD,  {Q1,Q1},  NLEV;
          FEATOOL,   QUAD,  {Q2,D1},  NLEV(1:end-1);
          FENICS,    TRI,   {P1,P1},  NLEV;
          FENICS,    TRI,   {P2,P1},  NLEV(1:end-1);
          OPENFOAM,  TRI,   {D0,D0},  NLEV;
          OPENFOAM,  QUAD,  {D0,D0},  NLEV;
          SU2,       TRI,   {P1,P1},  NLEV;
          SU2,       QUAD,  {Q1,Q1},  NLEV };


opt.basename = mfilename();
opt.solvers  = solvers;
opt.grids    = grids;
opt.cases    = cases;
opt.fcn_fea  = @l_setup_fea_struct;
opt.fcn_err  = @l_dragliftpres;
opt.fcn_proc = @l_process_data;

[ data, fea ] = run_featool_benchmarks( opt );


%------------------------------------------------------------------------------%
function [ fea ] = l_setup_fea_struct( grid_type, sfun, i_lev )

sf_u = sfun{1};
sf_p = sfun{2};

rho   = 1;           % Density.
miu   = 0.001;       % Molecular/dynamic viscosity.
umax  = 0.3;         % Maximum magnitude of inlet velocity.
umean = 2/3*umax;    % Mean inlet velocity.
h     = 0.41;        % Height of rectangular domain.
l     = 2.2;         % Length of rectangular domain.
xc    = 0.2;         % x-coordinate of cylinder center.
yc    = 0.2;         % y-coordinate of cylinder center.
diam  = 0.1;         % Diameter of cylinder.

% Geometry definition.
gobj1 = gobj_rectangle( 0, 2.2, 0, 0.41, 'R1' );
gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
fea.geom.objects = { gobj1 gobj2 };
fea = geom_apply_formula( fea, 'R1-C1' );
fea.sdim = { 'x', 'y' };

% Grid generation.
if( grid_type==1 || grid_type==2 )

  % Structured quadrilateral benchmark grid.
  ns = 8*2^(i_lev-1);
  r  = [0.05 0.06 0.08 0.11 0.15];
  x  = [0.41 0.5 0.7 1 1.4 1.8 2.2];
  for i=2:i_lev
    r = sort( [ r, (r(1:end-1)+r(2:end))/2 ] );
    x = sort( [ x, (x(1:end-1)+x(2:end))/2 ] );
  end

  grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
  grid2 = holegrid( ns, 2^(i_lev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
  grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
  grid3 = rectgrid( x, ns, [0.41 2.2;0 0.41] );
  fea.grid = gridmerge( grid3, 4, grid2, 6 );
  fea.grid.s(:) = 1;

  if( grid_type==2 )
    fea.grid = quad2tri( fea.grid );
  end
else
  fea.grid = gridgen( fea, 'hmax', abs(i_lev) );
end

% Boundary specifications.
DTOL      = sqrt(eps)*1e3;
i_inflow  = findbdr( fea, ['x<=',num2str(DTOL)] );     % Inflow boundary number.
i_outflow = findbdr( fea, ['x>=',num2str(l-DTOL)] );   % Outflow boundary number.
s_inflow  = ['4*',num2str(umax),'*(y*(',num2str(h),'-y))/',num2str(h),'^2'];   % Definition of inflow profile.
i_cyl     = findbdr( fea, ['sqrt((x-',num2str(xc),').^2+(y-',num2str(yc),').^2)<=(',num2str(diam/2+DTOL),')'] );    % Cylinder boundary number.

% Problem definition.
fea = addphys(fea,@navierstokes);
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
[fea.phys.ns.sfun{1:2}] = deal( sf_u );
fea.phys.ns.sfun{3} = sf_p;

% Boundary conditions.
fea.phys.ns.bdr.sel(i_inflow)  = 2;
fea.phys.ns.bdr.sel(i_outflow) = 4;
fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;

fea = parsephys(fea);
fea = parseprob(fea);

%------------------------------------------------------------------------------%
function [ err, cd, cl, dp ] = l_dragliftpres( fea )

I_CUB = 10;
rho = fea.phys.ns.eqn.coef{1,end}{1};
miu = fea.phys.ns.eqn.coef{2,end}{1};
umax  = 0.3;         % Maximum magnitude of inlet velocity.
umean = 2/3*umax;    % Mean inlet velocity.
xc    = 0.2;         % x-coordinate of cylinder center.
yc    = 0.2;         % y-coordinate of cylinder center.
diam  = 0.1;         % Diameter of cylinder.
DTOL  = sqrt(eps)*1e3;
i_cyl = findbdr( fea, ['sqrt((x-',num2str(xc),').^2+(y-',num2str(yc),').^2)<=(',num2str(diam/2+DTOL),')'] );    % Cylinder boundary number.

% Evaluate pressure difference.
dp = evalexpr('p',[0.15 0.25;0.2 0.2],fea);

% Calculate drag/lift (line integration method).
s_tfx = ['nx*p+',num2str(miu),'*(-2*nx*ux-ny*(uy+vx))'];
s_tfy = ['ny*p+',num2str(miu),'*(-nx*(vx+uy)-2*ny*vy)'];
s_cd  = ['2*(',s_tfx,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
s_cl  = ['2*(',s_tfy,')/(',num2str(rho),'*',num2str(umean),'^2*',num2str(diam),')'];
c_d1  = intbdr(s_cd,fea,i_cyl,I_CUB);
c_l1  = intbdr(s_cl,fea,i_cyl,I_CUB);


% Calculate drag/lift (volume integration method).
bdrm   = fea.bdr.bdrm{1};
ind_b  = [];
ind_bm = [];
for ii=i_cyl
  ind_b  = [ind_b, find(fea.grid.b(3,:)==ii)];
  ind_bm = [ind_bm, find(bdrm(3,:)==ii)];
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


fprintf('\n\nBenchmark quantities:\n\n')

fprintf('Drag coefficient,    cd = %6f (l), %6f (v) (Ref: 5.57953523384)\n',c_d1,c_d2)
fprintf('Lift coefficient,    cl = %6f (l), %6f (v) (Ref: 0.010618937712)\n',c_l1,c_l2)
fprintf('Pressure,            dp = %6f (Ref: 0.11752016697)\n',dp(1)-dp(2))

% Error computation.
cd  = [c_d1, c_d2];
cl  = [c_l1, c_l2];
dp  = dp(1) - dp(2);
err = [ abs(cd - 5.57953523384)/5.57953523384, ...
        abs(cl - 0.010618937712)/0.010618937712, ...
        abs(dp - 0.11752016697)/0.11752016697 ];
err = [ cd, cl, dp, err ];

%------------------------------------------------------------------------------%
function l_process_data( opt, data, fea )

COLORS  = {'r','g','b','c','m','y',[.7 .7 1],[.7 1 .7], ...
           [1 .7 .7],[1 0.5 0],[0.5 1 0],[1 0 0.5],'k'};
MARKERS = {'o','s','v','^','+','x','d','<','>','p','h','*'};

h_fig1 = figure();
h_fig2 = figure();
h_fig3 = figure();
for i=1:size(data,1)

  s_case = data{i,1};
  data_i = data{i,2};

  fprintf('\n\n%s\n|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n', s_case )
  fprintf('| ilev |   nel   |   nvt   |    ndof    |  t_sol  |  it |      cd_l      |      cd_v      |       cl_l      |       cl_v      |       dp       |\n' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
  fprintf('| %4i | %7i | %7i | %10i | %7.1f | %3i | %14.11f | %14.11f | %15.12f | %15.12f | %14.11f |\n', data_i(:,1:end-5).' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
  fprintf('| Ref. |         |         |            |         |     |  5.57953523384 |  5.57953523384 |  0.010618937712 |  0.010618937712 |  0.11752016697 |\n' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n\n\n' )

  figure(h_fig1)
  loglog( data_i(:,5), abs(data_i(:,end-4)), ...
          ['-',MARKERS{mod(i-1,length(MARKERS))+1}], ...
          'color', COLORS{mod(i-1,length(COLORS))+1}, ...
          'linewidth', 2 )
  hold on
  if( i==size(data,1) )
    legend( data(:,1), 'Location', 'best' )
    grid on
    xlabel( 'CPU time [s]' )
    ylabel( 'error(c_D)' )
    title(  'Cost vs. accuracy (drag coefficient)' )
  end

  figure(h_fig2)
  loglog( data_i(:,5), abs(data_i(:,end-2)), ...
          ['-',MARKERS{mod(i-1,length(MARKERS))+1}], ...
          'color', COLORS{mod(i-1,length(COLORS))+1}, ...
          'linewidth', 2 )
  hold on
  if( i==size(data,1) )
    legend( data(:,1), 'Location', 'best' )
    grid on
    xlabel( 'CPU time [s]' )
    ylabel( 'error(c_L)' )
    title(  'Cost vs. accuracy (lift coefficient)' )
  end

  figure(h_fig3)
  loglog( data_i(:,5), abs(data_i(:,end)), ...
          ['-',MARKERS{mod(i-1,length(MARKERS))+1}], ...
          'color', COLORS{mod(i-1,length(COLORS))+1}, ...
          'linewidth', 2 )
  hold on
  if( i==size(data,1) )
    legend( data(:,1), 'Location', 'best' )
    grid on
    xlabel( 'CPU time [s]' )
    ylabel( 'error(\Delta p)' )
    title(  'Cost vs. accuracy (pressure difference)' )
  end
end

figure(h_fig1)
print( '-r300', '-djpeg', fullfile(opt.workdir,'ns3_benchmark_drag') )
print( '-r300', '-dpng',  fullfile(opt.workdir,'ns3_benchmark_drag') )

figure(h_fig2)
print( '-r300', '-djpeg', fullfile(opt.workdir,'ns3_benchmark_lift') )
print( '-r300', '-dpng',  fullfile(opt.workdir,'ns3_benchmark_lift') )

figure(h_fig3)
print( '-r300', '-djpeg', fullfile(opt.workdir,'ns3_benchmark_pres') )
print( '-r300', '-dpng',  fullfile(opt.workdir,'ns3_benchmark_pres') )
