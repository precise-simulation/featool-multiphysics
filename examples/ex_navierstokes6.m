function [ fea, out ] = ex_navierstokes6( varargin )
%EX_NAVIERSTOKES6 2D incompressible time-dependent flow around a cylinder.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES6( VARARGIN ) Benchmark example for time-dependent flow around a cylinder.
%
%   References:
%
%   [1] John V. Reference values for drag and lift of a two-dimensional
%       time-dependent flow around a cylinder. International Journal for
%       Numerical Methods in Fluids 2004; 44:777-788.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {0.001}         Molecular/dynamic viscosity
%       igrid       scalar {2}             Grid type: >0 regular (igrid refinements)
%                                                     <0 unstruc. grid (with hmax=|igrid|)
%       dt          scalar {0.02}          Time step size
%       ischeme     scalar {3}             Time stepping scheme
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sf_disc1}      Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2019 Precise Simulation, Ltd.


cOptDef = { ...
  'rho',      1; ...
  'miu',      0.001; ...
  'igrid',    2; ...
  'dt',       0.02; ...
  'ischeme'   3; ...
  'sf_u',     'sflag2'; ...
  'sf_p',     'sf_disc1'; ...
  'iplot',    1; ...
  'itest',    0; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Grid generation.
fea.sdim = { 'x' 'y' };
if( opt.igrid>=1 )
  fea.grid = cylbenchgrid( opt.igrid );
else
  gobj1 = gobj_rectangle( 0, 2.2, 0, 0.41, 'R1' );
  gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
  fea.geom.objects = { gobj1 gobj2 };
  fea = geom_apply_formula( fea, 'R1-C1' );
  fea.grid = gridgen( fea, 'hmax', abs(opt.igrid), 'fid', fid );
end


% Boundary conditions.
inflow_boundary  = findbdr( fea, ['x<=',num2str(sqrt(eps))] );   % Inflow boundary number.
outflow_boudnary = findbdr( fea, ['x>=',num2str(2.1)] );         % Outflow boundary number.
inflow_bc = '6*sin(pi*t/8)*(y*(0.41-y))/0.41^2';


% Problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun = { opt.sf_u opt.sf_u opt.sf_p };
fea.phys.ns.bdr.sel( inflow_boundary  ) = 2;
fea.phys.ns.bdr.sel( outflow_boudnary ) = 3;
fea.phys.ns.bdr.coef{2,end}{1,inflow_boundary} = inflow_bc;
fea = parsephys( fea );


% Call solver
[fea.sol.u,tlist] = solvetime( fea, 'fid', fid, 'tstep', opt.dt, 'tmax', 8*(~opt.itest), 'maxnit', 1+4*(~opt.itest), 'ischeme', opt.ischeme );


% Benchmark quantities.
[ c_d, c_l, dp ] = calc_bench_quant( fea, 1 );
[c_d_max, i] = max( c_d );
t_c_d_max    = tlist( i );
[c_l_max, i] = max( c_l );
t_c_l_max    = tlist( i );
[~, i]       = min(abs(tlist-8));
dp_t8        = dp(i);
out.c_d = c_d;
out.c_l = c_l;
out.dp  = dp;
out.c_d_max = c_d_max;
out.c_l_max = c_l_max;
out.t_c_d_max = t_c_d_max;
out.t_c_l_max = t_c_l_max;
out.dp_t8 = dp_t8;

if( ~isempty(fid) )
  fmtf = '%12.8f |';
  fmts = '%12s |';
  fmt  = ['|      ',repmat(fmtf,1,5),'\n'];
  fmts = ['|      ',repmat(fmts,1,5),'\n'];
  fmtl = ['|------',repmat('-------------+',1,5)];
  fmtl = [fmtl(1:end-1),'|\n'];
  fprintf( fid, '\n\n' );
  fprintf( fid, fmtl );
  fprintf( fid, fmts, 't(cd_max)', 'cd_max', 't(cl_max)', 'cl_max', 'dp(t=8)' );
  fprintf( fid, fmtl );
  fprintf( fid, fmt, t_c_d_max, c_d_max, t_c_l_max, c_l_max, dp_t8 );
  fprintf( fid, fmtl );
  ref_data = [ 3.93625 2.950923849 5.6925 0.47834818 -0.11162153 ];
  fmt  = ['| Ref. ',repmat(fmtf,1,5),'\n'];
  fprintf( fid, fmt, ref_data );
  fprintf( fid, fmtl );
  fprintf( fid, '\n\n' );
end


% Postprocessing.
if( opt.iplot>0 )
  figure
  fea = parseprob( fea );

  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'}, 'arrowcolor', 'k' )
  title('Velocity field at t=8')

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


if( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ c_d, c_l, dp ] = calc_bench_quant( fea, umean )
% CALC_BENCH_QUANT Calculate benchmark quantities

  i_type  = 2;   % 1-line integration, 2-volume integration
  i_cub   = 10;  % Quadrature order for line integration.
  i_alpha = 1;   % Shape function orde for alpha 'a' field.
  cylinder_boundary = findbdr( fea, 'sqrt((x-0.2).^2+(y-0.2).^2)<0.1' );   % Cylinder boundary number.
  rho = fea.coef{1,end}{1};
  miu = fea.coef{2,end}{1};

  if( i_type==1 )

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
    if( i_alpha==1 )
      fea.dvar = [ fea.dvar {'a'}      ];
      fea.sfun = [ fea.sfun {'sflag1'} ];
      fea      = parseprob(fea);
      u_a      = zeros(size(fea.grid.p,2),1);
      ind_b    = find(ismember(fea.grid.b(3,:),cylinder_boundary));
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
      for ii=cylinder_boundary
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
  for i_sol=1:size(u,2);
    fea.sol.u = u(:,i_sol);

    if( i_type==1 )   % Line integration.
      c_d(i_sol) = intbdr( s_cd_line, fea, cylinder_boundary, i_cub );
      c_l(i_sol) = intbdr( s_cl_line, fea, cylinder_boundary, i_cub );
    else   % Volume integration.
      fea.sol.u = [fea.sol.u;u_a];
      fea.eqn   = struct;
      fea.bdr   = struct;
      fea       = parseprob(fea);

      c_d(i_sol) = intsubd( s_cd_vol, fea, [], [], 3 );
      c_l(i_sol) = intsubd( s_cl_vol, fea, [], [], 3 );
    end

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
