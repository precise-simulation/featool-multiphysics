function [ fea, out ] = ex_navierstokes6( varargin )
%EX_NAVIERSTOKES6 Instationary laminar 2D flow around a cylinder (Re(t)<=100).
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES3( VARARGIN ) 2D validation and CFD
%   benchmark for instationary laminar incompressible flow around a
%   cylinder (Reynolds number, Re(t) <= 100).
%
%   The maximum drag and lift coefficients (with corresponding
%   incidence times), and pressure difference between the front and
%   back of the cylinder at t = 8s are computed and compared with
%   reference values [1].
%
%   Reference:
%
%   [1] John V. Reference values for drag and lift of a two-dimensional
%       time-dependent flow around a cylinder. International Journal for
%       Numerical Methods in Fluids, 44:777-788, 2004.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar {2}             Grid type: >0 regular (igrid refinements)
%                                                     <0 unstruc. grid (with hmax=|igrid|)
%       dt          scalar {0.02}          Time step size
%       instatbc    boolean {true}         Use instationary boundary conditions
%       ischeme     scalar {2}             Time stepping scheme
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sf_disc1}      Shape function for pressure
%       solver      string {default}       Solver selection default, openfoam, or su2
%       iplot       logical false/{true}   Plot and visualize solution
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES3, EX_NAVIERSTOKES13

% Copyright 2013-2026 Precise Simulation, Ltd.


cOptDef = { 'igrid',    2;
            'dt',       0.02;
            'instatbc'  true;
            'ischeme'   2;
            'sf_u',     'sflag2';
            'sf_p',     'sf_disc1';
            'solver',   '';
            'iplot',    true;
            'tol',      [0.01, 0.01, 0.1, 0.5, 0.2];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
if ( any(strcmpi(opt.solver,{'SU2'})) )
  if ( ~isempty(opt.fid) )
    warning('SU2 does not support intationary BCs.')
  end
  opt.instatbc = false;
end
if ( strcmpi(opt.solver,'OPENFOAM') && ~got.igrid )
  opt.igrid = 3;
end

% Model parameters.
rho = 1;      % Density.
miu = 0.001;  % Molecular/dynamic viscosity.

if ( opt.instatbc )
             % t_c_d_max c_d_max t_c_l_max c_l_max dp_t8
  ref_data = [3.93625, 2.950923849, 5.6925, 0.47834818, -0.11162153];
else
             % c_d_max c_l_max St dp_f2
  ref_data = [3.23, 1.00, 0.3, 2.48];
  if( ~got.tol )
    opt.tol = 1;
  end
end

% Geometry.
h    = 0.41;  % Height of rectangular domain.
l    = 2.2;   % Length of rectangular domain.
xc   = 0.2;   % x-coordinate of cylinder center.
yc   = 0.2;   % y-coordinate of cylinder center.
diam = 0.1;   % Diameter of cylinder.

fea.sdim = { 'x', 'y' };
R1 = gobj_rectangle( 0, l, 0, h, 'R1' );
C1 = gobj_circle( [xc, yc], diam/2, 'C1' );
fea.geom.objects = { R1, C1 };
fea = geom_apply_formula( fea, 'R1-C1' );


% Grid/mesh generation.
if ( opt.igrid >= 1 )
  fea.grid = l_cylbenchgrid_2d( opt.igrid );
else
  fea.grid = gridgen( fea, 'hmax', abs(opt.igrid), 'fid', opt.fid );
end


% Boundary identification.
dtol = sqrt(eps) * 1e3;
ind_inflow = findbdr(fea, sprintf('x <= %.16g', dtol));  % Inflow boundary number.
ind_outflow = findbdr(fea, sprintf('x >= (%.16g - %.16g)', l, dtol ));  % Outflow boundary number.
ind_cylinder = findbdr(fea, sprintf('sqrt((x-%.16g).^2+(y-%.16g).^2) <= %.16g', xc, yc, diam/2 + dtol));  % Cylinder boundary numbers.


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
if ( any(strcmpi(opt.solver,{'OPENFOAM','SU2'})) )
  fea.phys.ns.sfun = { 'sflag1', 'sflag1', 'sflag1' };
else
  fea.phys.ns.sfun = { opt.sf_u, opt.sf_u, opt.sf_p };
end

% Boundary condition definitions.
fea.phys.ns.bdr.sel(ind_inflow) = 2;
fea.phys.ns.bdr.sel(ind_outflow) = 4;
if ( opt.instatbc )
  s_inflow = '6*sin(pi*t/8)*(y*(0.41-y))/0.41^2';
else
  s_inflow = '6*(y*(0.41-y))/0.41^2';
end
fea.phys.ns.bdr.coef{2,end}{1,ind_inflow} = s_inflow;


% Call solver.
fea = parsephys(fea);
fea = parseprob(fea);  % Check and parse problem struct.
switch ( upper(opt.solver) )

  case 'OPENFOAM'
    if ( opt.ischeme == 1 )
      ddtSchemes = 'backward';
    else
      ddtSchemes = 'CrankNicolson 0.9';
    end
    fid = opt.fid; if ( ~got.fid ), fid = []; end
    [fea.sol.u,tlist] = openfoam( fea, 'application', 'pimpleFoam', 'maxCo', 1.0, 'ddtSchemes', ddtSchemes, ...
                                  'deltaT', opt.dt, 'endTime', 8, 'nproc', 1, 'fid', fid, 'logfid', opt.fid );
    fea.sol.u = fea.sol.u(:,tlist > 0);
    tlist = tlist(tlist > 0);

  case 'SU2'
    fid = opt.fid; if ( ~got.fid ), fid = []; end
    [fea.sol.u,tlist] = su2( fea, 'fid', fid, 'logfid', opt.fid, 'tstep', opt.dt, 'tmax', 8, 'ischeme', opt.ischeme, 'wrtfreq', floor(8/opt.dt/400) );

  otherwise
    [fea.sol.u,tlist] = solvetime( fea, 'fid', opt.fid, 'tstep', opt.dt, 'tmax', 8, 'maxnit', 5, 'ischeme', opt.ischeme );
end

% Benchmark quantities.
[ c_d, c_l, dp ] = calc_bench_quants( fea, 1, ind_cylinder );
[c_d_max, i] = max( c_d );
t_c_d_max    = tlist( i );
[c_l_max, i] = max( c_l );
t_c_l_max    = tlist( i );
out.c_d = c_d;
out.c_l = c_l;
out.dp  = dp;
out.c_d_max = c_d_max;
out.c_l_max = c_l_max;
out.t_c_d_max = t_c_d_max;
out.t_c_l_max = t_c_l_max;
if ( opt.instatbc )
  [~, i] = min(abs(tlist-8));
  dp_t8  = dp(i);
  out.dp_t8 = dp_t8;

  comp_data = [ t_c_d_max c_d_max t_c_l_max c_l_max dp_t8 ];
  out.err = abs(ref_data-comp_data)./abs(ref_data);
  out.pass = all( out.err < opt.tol );

else
  found = false;
  for j=i+1:length(c_l)
    if ( j==length(c_l) )
      break
    end
    if ( j>1 && c_l(j) > c_l(j-1) && c_l(j) > c_l(j+1) )
      found = true;
      c_l_max2 = c_l(j);
      t_c_l_max2 = tlist(j);
      break
    end
  end
  for j=i-1:-1:1
    if ( j>1 && c_l(j) > c_l(j-1) && c_l(j) > c_l(j+1) )
      found = true;
      c_l_max2 = c_l(j);
      t_c_l_max2 = tlist(j);
      break
    end
  end
  f = 1/(t_c_l_max2-t_c_l_max);
  St = 0.1*f/1.5;

  out.St = St;
  t_dp = t_c_l_max + 1/2/f;
  [~,j] = find( tlist > t_dp, 1 );
  dp_f2 = dp(j-1) + (t_dp - tlist(j-1))/(tlist(j) - tlist(j-1)) * (dp(j)-dp(j-1));
  out.dp_f2 = dp_f2;

  c_d(tlist < 0.5) = 0;
  [c_d_max, i] = max( c_d );
  comp_data = [ c_d_max c_l_max St dp_f2 ];
  out.err = abs(ref_data-comp_data)./abs(ref_data);
  out.pass = all( out.err([1,2,4]) < opt.tol );
end


if ( ~isempty(opt.fid) )
  fmtf = '%12.8f |';
  fmts = '%12s |';
  fmt  = ['|      ',repmat(fmtf,1,4+double(opt.instatbc)),'\n'];
  fmts = ['|      ',repmat(fmts,1,4+double(opt.instatbc)),'\n'];
  fmtl = ['|------',repmat('-------------+',1,4+double(opt.instatbc))];
  fmtl = [fmtl(1:end-1),'|\n'];
  fprintf( opt.fid, '\n\n' );
  fprintf( opt.fid, fmtl );
  if ( opt.instatbc )
    fprintf( opt.fid, fmts, 't(cd_max)', 'cd_max', 't(cl_max)', 'cl_max', 'dp(t=8)' );
    fprintf( opt.fid, fmtl );
    fprintf( opt.fid, fmt, t_c_d_max, c_d_max, t_c_l_max, c_l_max, dp_t8 );
  else
    fprintf( opt.fid, fmts, 'cd_max', 'cl_max', 'St', 'dp' );
    fprintf( opt.fid, fmtl );
    fprintf( opt.fid, fmt, c_d_max, c_l_max, St, dp_f2 );
  end
  fprintf( opt.fid, fmtl );
  fmt = ['| Ref. ',repmat(fmtf,1,4+double(opt.instatbc)),'\n'];
  fprintf( opt.fid, fmt, ref_data );
  fprintf( opt.fid, fmtl );
  fprintf( opt.fid, '\n\n' );
end


% Postprocessing.
if ( opt.iplot )
  figure
  fea = parseprob( fea );

  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'}, 'arrowcolor', 'k' )
  title('Velocity field at t=8')

  ix = tlist > 1;
  if ( any(ix) )
    tlist = tlist(ix);
    c_d = c_d(ix);
    c_l = c_l(ix);
    dp  = dp(ix);
  end

  subplot(2,2,3)
  plot( tlist, c_d )
  title('drag coefficient')

  subplot(2,2,4)
  plot( tlist, c_l )
  title('lift coefficient')

  subplot(2,2,2)
  plot( tlist, dp )
  title('pressure difference')
end


if ( ~nargout )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ c_d, c_l, dp ] = calc_bench_quants( fea, umean, ind_cylinder )
% CALC_BENCH_QUANT Calculate benchmark quantities

i_type  = 2;   % 1-line integration, 2-volume integration
i_cub   = 10;  % Quadrature order for line integration.
i_alpha = 1;   % Shape function orde for alpha 'a' field.
rho = fea.coef{1,end}{1};
miu = fea.coef{2,end}{1};

if ( i_type == 1 )

  s_tfx_line = ['nx*p+',num2str(miu),'*(-2*nx*ux-ny*(uy+vx))'];
  s_tfy_line = ['ny*p+',num2str(miu),'*(-nx*(vx+uy)-2*ny*vy)'];
  s_cd_line  = ['2*(',s_tfx_line,')/(',num2str(rho),'*',num2str(umean),'^2*0.1)'];
  s_cl_line  = ['2*(',s_tfy_line,')/(',num2str(rho),'*',num2str(umean),'^2*0.1)'];

else

  s_tfx_vol  = ['ax*p+',num2str(miu),'*(-2*ax*ux-ay*(uy+vx))-(u*ux+v*uy)*a'];
  s_tfy_vol  = ['ay*p+',num2str(miu),'*(-ax*(vx+uy)-2*ay*vy)-(u*vx+v*vy)*a'];
  s_cd_vol   = ['2*(',s_tfx_vol,')/(',num2str(rho),'*',num2str(umean),'^2*0.1)'];
  s_cl_vol   = ['2*(',s_tfy_vol,')/(',num2str(rho),'*',num2str(umean),'^2*0.1)'];

  % Extend fea struct with an alpha field 'a' with values one on the cylinder
  % and zero everywhere else (for volume integraion).
  if ( i_alpha == 1 )
    fea.dvar = [ fea.dvar {'a'}      ];
    fea.sfun = [ fea.sfun {'sflag1'} ];
    fea      = parseprob(fea);
    u_a      = zeros(size(fea.grid.p,2),1);
    ind_b    = find(ismember(fea.grid.b(3,:),ind_cylinder));
    ind_c    = fea.grid.b(1,ind_b);
    ind_ei   = fea.grid.b(2,ind_b);
    ind_ej   = mod(ind_ei,size(fea.grid.c,1)) + 1;
    ind_c    = sub2ind( size(fea.grid.c), [ind_ei;ind_ej], [ind_c;ind_c] );
    ind_v    = fea.grid.c(ind_c);
    u_a(ind_v) = 1;
  else
    bdrm   = fea.bdr.bdrm{1};
    ind_b  = [];
    ind_bm = [];
    for ii=ind_cylinder
      ind_b  = [ind_b  find(fea.grid.b(3,:)==ii)];
      ind_bm = [ind_bm find(bdrm(3,:)==ii)];
    end
    ind_c    = fea.grid.b(1,ind_b);
    ind_gdof = bdrm(4,ind_bm);

    fea.dvar = [ fea.dvar {'a'}      ];
    fea.sfun = [ fea.sfun {'sflag2'} ];
    fea      = parseprob(fea);
    n_dof    = max(fea.eqn.dofm{1}(:));
    u_a      = zeros(n_dof,1);
    u_a(ind_gdof) = 1;
  end

end

% Loop over all time steps and compute.
u = fea.sol.u;
if ( i_type == 1 )   % Line integration.
  c_expr = {s_cd_line, s_cl_line};
else
  c_expr = {s_cd_vol, s_cl_vol};
end

c_d = zeros(1,size(u,2));
c_l = zeros(1,size(u,2));
for i=1:length(c_expr)
  for i_sol=1:size(u,2);
    fea.sol.u = u(:,i_sol);

    if ( i_type == 1 )   % Line integration.
      val = intbdr( c_expr{i}, fea, ind_cylinder, i_cub );

    else   % Volume integration.
      fea.sol.u = [fea.sol.u;u_a];
      fea.eqn = struct; fea.bdr = struct;
      fea = parseprob(fea);

      val = intsubd( c_expr{i}, fea, [], [], 3 );
    end

    if ( i == 1 )
      c_d(i_sol) = val;
    else
      c_l(i_sol) = val;
    end
  end
end

for i_sol=1:size(u,2);
  fea.sol.u = u(:,i_sol);

  % Pressure difference.
  p_i = [nan nan];
  x = 0.15;
  while( isnan(p_i(1)) )
    p_i(1) = evalexpr('p',[x;0.2],fea);
    x = x - 0.005;
  end
  x = 0.25;
  while( isnan(p_i(2)) )
    p_i(2) = evalexpr('p',[x;0.2],fea);
    x = x + 0.005;
  end
  dp(i_sol) = p_i(1) - p_i(2);

end

%------------------------------------------------------------------------------%
function [ grid ] = l_cylbenchgrid_2d( nlev )
% CYLBENCHGRID_2D Generate quadrilateral grid for the 2D DFG cylinder benchmark.

ns = 8 * 2^(nlev-1);
r  = [0.05, 0.06, 0.08, 0.11, 0.15];
x  = [0.41, 0.5, 0.7, 1, 1.4, 1.8, 2.2];
for ilev=2:nlev
  r = sort( [ r (r(1:end-1)+r(2:end))/2 ] );
  x = sort( [ x (x(1:end-1)+x(2:end))/2 ] );
end

grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
grid2 = holegrid( ns, 2^(nlev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
grid3 = rectgrid( x, ns, [0.41, 2.2; 0, 0.41] );
grid  = gridmerge( grid3, 4, grid2, 6 );

grid.s(:) = 1;
