function [ data, fea ] = ex_navierstokes3b( varargin )
%EX_NAVIERSTOKES3B 2D Example for incompressible stationary flow around a cylinder.
%
%   [ DATA ] = EX_NAVIERSTOKES3B Flow around a cylinder CFD benchmark
%
%   CFD benchmarking script comparing the OpenFOAM CFD solver to the
%   fully coupled and monolithic FEniCS and FEATool Multiphysics solvers,
%   for a 2D stationary and incompressible flow around a cylinder test case.
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
%   See also EX_NAVIERSTOKES3.

% Copyright 2013-2019 Precise Simulation, Ltd.


% Simulation parameters.
nlev   = 1:3;        % Grid levels.
nsolv  = 1:3;        % Solvers 1: FEATool, 2: FEniCS, 3: OpenFOAM
iplot  = 1;          % Plot and visualize results.
isave  = true;       % Save benchmark results to file.
iimg   = false;      % Save images.
ldir   = fullfile(tempdir(),'ns3_bench');   % Benchmark directory.
lfile  = 'ns3_bench.log';   % Logfile name.

% Test case settings.
rho       = 1;           % Density.
miu       = 0.001;       % Molecular/dynamic viscosity.
umax      = 0.3;         % Maximum magnitude of inlet velocity.
umean     = 2/3*umax;    % Mean inlet velocity.
% Geometry.
h         = 0.41;        % Height of rectangular domain.
l         = 2.2;         % Length of rectangular domain.
xc        = 0.2;         % x-coordinate of cylinder center.
yc        = 0.2;         % y-coordinate of cylinder center.
diam      = 0.1;         % Diameter of cylinder.
% FEM discretization parameters.
sf_u      = 'sflag2';    % FEM shape function type for velocity.
sf_p      = 'sflag1';    % FEM shape function type for pressure.


% Run benchmarks.
rmkdir(ldir);
solvers = {'FEATool','FEniCS','OpenFOAM'};
data = {};
for isolv=nsolv
  data_i = [];
  for ilev=nlev
    clear fea

    % Geometry definition.
    gobj1 = gobj_rectangle( 0, 2.2, 0, 0.41, 'R1' );
    gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
    fea.geom.objects = { gobj1 gobj2 };
    fea = geom_apply_formula( fea, 'R1-C1' );
    fea.sdim = { 'x', 'y' };

    % Grid generation.
    if( ilev>=1 )

      % Structured quadrilateral benchmark grid.
      ns = 8*2^(ilev-1);
      r  = [0.05 0.06 0.08 0.11 0.15];
      x  = [0.41 0.5 0.7 1 1.4 1.8 2.2];
      for i=2:ilev
        r = sort( [ r, (r(1:end-1)+r(2:end))/2 ] );
        x = sort( [ x, (x(1:end-1)+x(2:end))/2 ] );
      end

      grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
      grid2 = holegrid( ns, 2^(ilev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
      grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
      grid3 = rectgrid( x, ns, [0.41 2.2;0 0.41] );
      fea.grid = gridmerge( grid3, 4, grid2, 6 );
      fea.grid.s(:) = 1;

      if( isolv==2 && size(fea.grid.c,1)==4 )
        fea.grid = quad2tri( fea.grid );
      end
    else
      fea.grid = gridgen( fea, 'hmax', abs(ilev) );
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
    fea.phys.ns.sfun            = { sf_u, sf_u, sf_p };

    % Boundary conditions.
    fea.phys.ns.bdr.sel(i_inflow)  = 2;
    fea.phys.ns.bdr.sel(i_outflow) = 4;
    fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;

    fprintf( 1, '\n%s - Level %i : %g\n\n', solvers{isolv}, ilev );

    % Parse and solve problem.
    fea = parsephys(fea);
    if( isolv==2 )
      fea.dvar{3} = 'sflag1';
    end
    fea = parseprob(fea);
    fid = [];
    switch( isolv )

      case 1   % FEATool (analytic Newton Jacobian)
        jac.form  = { [1;1] [1;1] [] ;
                      [1;1] [1;1] [] ;
                      []    []    [] };
        jac.coef  = { [num2str(rho),'*ux'] [num2str(rho),'*uy'] [] ;
                      [num2str(rho),'*vx'] [num2str(rho),'*vy'] [] ;
                      []                   []                   [] };

        fid = fopen( fullfile(ldir,lfile), 'w+' );
        fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac );

      case 2   % Use external FEniCS solver.
        fea = fenics( fea, 'fdir', ldir, 'fname', lfile, 'clear', false );

      case 3   % Use external OpenFOAM CFD solver.
        fea.sol.u = openfoam( fea, 'casedir', ldir, 'logfname', lfile, 'clean', false );
    end
    [t_solv,it] = l_parse_logfile( isolv, ldir, lfile, fid );


    % Calculate benchmark quantities and error.
    [cd,cl,dp,err] = l_dragliftpres( fea, rho, miu, umean, diam, i_cyl );

    nel  = size(fea.grid.c,2);
    nvt  = size(fea.grid.p,2);
    ndof = sum(fea.eqn.ndof);
    if( isolv==3 )
      ndof = 3*nel;   % OpenFOAM uses cell centered dofs.
    end
    data_i = [ data_i; [ ilev, nel, nvt, ndof, t_solv, it, cd, cl, dp, err ] ];

    delete(fullfile(ldir,lfile))
  end

  data = [ data; { data_i } ];
  if( isave )
    save 'ns3_benchmark_data' data
  end
end


% Data processing.
solvers  = solvers(nsolv);
marker   = {'o','s','v'};
for i=1:length(data)

  data_i = data{i};

  fprintf('\n\n%s\n|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n', solvers{i} )
  fprintf('| ilev |   nel   |   nvt   |    ndof    |  t_sol  |  it |      cd_l      |      cd_v      |       cl_l      |       cl_v      |       dp       |\n' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
  fprintf('| %4i | %7i | %7i | %10i | %7.1f | %3i | %14.11f | %14.11f | %15.12f | %15.12f | %14.11f |\n', data_i(:,1:end-5).' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
  fprintf('| Ref. |         |         |            |         |     |  5.57953523384 |  5.57953523384 |  0.010618937712 |  0.010618937712 |  0.11752016697 |\n' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n\n\n' )

  if( iplot )
    figure(1)
    loglog( data_i(:,5), data_i(:,end-4), ['-',marker{i}], 'linewidth', 2 )
    hold on
    if( i==length(data) )
      legend( solvers{:}, 'Location', 'SouthWest' )
      grid on
      xlabel( 'CPU time [s]' )
      ylabel( 'error(c_D)' )
      title( 'Cost vs. accuracy (drag coefficient)' )
    end

    figure(2)
    loglog( data_i(:,5), data_i(:,end-2), ['-',marker{i}], 'linewidth', 2 )
    hold on
    if( i==length(data) )
      legend( solvers{:} )
      grid on
      xlabel( 'CPU time [s]' )
      ylabel( 'error(c_L)' )
      title( 'Cost vs. accuracy (lift coefficient)' )
    end

    figure(3)
    loglog( data_i(:,5), data_i(:,end), ['-',marker{i}], 'linewidth', 2 )
    hold on
    if( i==length(data) )
      legend( solvers{:} )
      grid on
      xlabel( 'CPU time [s]' )
      ylabel( 'error(\Delta p)' )
      title( 'Cost vs. accuracy (pressure difference)' )
    end
  end
end
if( iplot && iimg )
  figure(1)
  print -r300 -dpng ns3_benchmark_drag
  figure(2)
  print -r300 -dpng ns3_benchmark_lift
  figure(3)
  print -r300 -dpng ns3_benchmark_pres
end


% Flow field visualization.
if( iplot>1 )
  figure
  subplot(2,1,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u','v'} )
  title( 'Velocity field' )
  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'p' )
  title( 'Pressure' )

  if( iimg )
    print -r300 -dpng ns3_benchmark_vis
  end
end

if( ~nargout )
  clear data fea
end



%------------------------------------------------------------------------------%
function [ cd, cl, dp, err ] = l_dragliftpres( fea, rho, miu, umean, diam, i_cyl )

I_CUB = 10;

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

%------------------------------------------------------------------------------%
function [t_solv,it] = l_parse_logfile( isolv, ldir, lfile, fid )
% Parse log file to determine solver time and number of iterations.

switch( isolv )

  case 1   % FEATool
    frewind( fid );
    sline = '';
    while( ~any(findstr(sline,'it')) )
      sline = fgetl( fid );
    end
    sline = fgetl( fid );
    sline = fgetl( fid );
    it = 0;
    while 1
      sline = fgetl( fid );
      if( any(findstr(sline,'---')) )
        break
      end
      it = it + 1;
    end
    while( ~any(findstr(sline,'t_tot')) )
      sline = fgetl( fid );
    end
    t_solv = str2num(sline(13:end));

  case 2   % FEniCS
    fid = fopen( fullfile(ldir,lfile), 'r' );
    while 1
      sline = fgetl( fid );
      if ~ischar(sline), break, end
      if( any(findstr(sline,'t_solve')) )
        t_solv = str2num(sline(10:end));
      end
      if( any(findstr(sline,'Newton iteration')) )
        ix = find(sline == ':');
        c = regexp( strtrim(sline(1:(ix(1)-1))), '\s+', 'split' );
        for i=1:length(c)
          it = str2num(c{i});
          if( ~isempty(it) )
            break
          end
        end
      end

    end

  case 3   % OpenFOAM
    fid = fopen( fullfile(ldir,lfile), 'r' );
    while 1
      sline = fgetl( fid );
      if ~ischar(sline), break, end
      if( any(findstr(sline,'ExecutionTime')) )
        c = regexp( strtrim(sline(17:end)), '\s+', 'split' );
        t_solv = str2num(c{1});
      end
      if( length(sline)>7 && strcmp(sline(1:7),'Time = ') )
        c = regexp( strtrim(sline), '\s+', 'split' );
        for i=1:length(c)
          it = str2num(c{i});
          if( ~isempty(it) )
            break
          end
        end
      end
    end

end
fclose( fid );
