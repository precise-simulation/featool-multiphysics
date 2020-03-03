function [ data, fea ] = ex_navierstokes8b( varargin )
%EX_NAVIERSTOKES8B Axisymmtric laminar pipe flow benchmark.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES8B Sets up and solves stationary
%   axisymmetric laminar Poiseuille flow in a circular pipe. The
%   inflow profile is constant and the outflow should assume a
%   parabolic profile. Compares the FEATool, FEniCS, and OpenFOAM
%   solvers.
%
%   See also EX_NAVIERSTOKES8.

% Copyright 2013-2020 Precise Simulation, Ltd.


% Simulation parameters.
nlev   = 1:3;        % Grid levels.
nsolv  = 1:3;        % Solvers 1: FEATool, 2: FEniCS, 3: OpenFOAM
igrid  = 1;          % 1: Structured grid, 0: Unstructured grid
iplot  = true;       % Plot and visualize results.
isave  = true;       % Save benchmark results to file.
iimg   = false;      % Save images.
ldir   = fullfile(tempdir(),'ns8_bench');   % Benchmark directory.
lfile  = 'ns8_bench.log';   % Logfile name.

% Test case settings.
rho    = 2;          % Density.
miu    = 3;          % Viscosity.
w_in   = 1;          % Inlet velocity.
% Geometry and grid parameters.
r      = 1;          % Pipe radius.
l      = 5;          % Pipe length.
hmax   = 1/5;        % Grid size on level 1.
% FEM discretization parameters.
isobdr = true;       % Enforce stright outflow boundary condition.
sf_u   = 'sflag2';   % FEM shape function for velocity.
sf_p   = 'sflag1';   % FEM shape function for velocity.


% Run benchmarks.
rmkdir(ldir);
solvers = {'FEATool','FEniCS','OpenFOAM'};
data = {};
for isolv=nsolv
  data_i = [];
  vel_profile = cell(0,4);
  for ilev=nlev
    clear fea

    % Geometry definition.
    gobj1 = gobj_rectangle( 0, r, 0, l, 'R1' );
    fea.geom.objects = { gobj1 };
    fea.sdim = { 'r' 'z' };

    % Grid generation.
    hmax_i = hmax/2^(ilev-1);
    if( igrid==1 )
      fea.grid = quad2tri(rectgrid(round(r/hmax_i),round(l/hmax_i),[0 r;0 l]));
    else
      fea.grid = gridgen( fea, 'hmax', hmax_i );
    end

    % Boundary specifications.
    i_inflow   = 1;   % Inflow boundary number.
    i_outflow  = 3;   % Outflow boundary number.
    i_symmetry = 4;   % Symmetry axis boundary number.

    % Problem definition.
    IS_AXI = true;
    fea = addphys( fea, {@navierstokes, IS_AXI} );
    fea.phys.ns.eqn.coef{1,end} = { rho };
    fea.phys.ns.eqn.coef{2,end} = { miu };
    fea.phys.ns.sfun            = { sf_u sf_u sf_p };

    % Boundary conditions.
    fea.phys.ns.bdr.sel(i_inflow)   = 2;
    fea.phys.ns.bdr.sel(i_outflow)  = 4;
    fea.phys.ns.bdr.coef{2,end}{2,i_inflow} = w_in;
    fea.phys.ns.bdr.sel(i_symmetry) = 5;

    % Parse problem.
    fea = parsephys(fea);
    fea = parseprob(fea);
    if( isobdr )
      % Enforce straight-outflow (u=0).
      fea.bdr.d{1}{i_outflow} = 0;
      fea.bdr.d{2}{i_outflow} = [];
      fea.bdr.n{1}{i_outflow} = [];
      fea.bdr.n{2}{i_outflow} = 0;

      % Enforce slip on symmetry axis (u=0).
      fea.bdr.d{1}{i_symmetry} = 0;
      fea.bdr.d{2}{i_symmetry} = [];
      fea.bdr.n{1}{i_symmetry} = [];
      fea.bdr.n{2}{i_symmetry} = 0;
    end

    % Solve problem.
    fid = [];
    switch( isolv )

      case 1   % FEATool (analytic Newton Jacobian)
        jac.form  = {[1;1] [1;1] [];
                     [1;1] [1;1] [];
                     []    []    []};
        jac.coef  = {'r*rho_ns*ur' 'r*rho_ns*uz' [];
                     'r*rho_ns*wr' 'r*rho_ns*wz' [];
                     []            []            []};
        fid = fopen( fullfile(ldir,lfile), 'w+' );
        fea.sol.u = solvestat( fea, 'nsolve', 2, 'jac', jac, 'fid', fid );

      case 2   % FEniCS
        fea = fenics( fea, 'fdir', ldir, 'fname', lfile, 'clear', false );

      case 3   % OpenFOAM
        fea.sol.u = openfoam( fea, 'casedir', ldir, 'logfname', lfile, 'clean', false );
    end
    [t_solv,it] = l_parse_logfile( isolv, ldir, lfile, fid );


    % Error checking.
    N_P = 20;
    p_r   = linspace( 0, r, N_P );
    p_z   = 0.9*l*ones( 1, N_P );
    w_vel = evalexpr( 'w', [p_r;p_z], fea )';
    w_ref = 2*w_in*( 1 - (p_r/r).^2 );   % Analytic solution 2*V_avg*((1-(r/R)^2).
    err   = sqrt( sum((w_vel-w_ref).^2)/sum(w_ref.^2) );

    nel  = size(fea.grid.c,2);
    nvt  = size(fea.grid.p,1);
    ndof = sum(fea.eqn.ndof);
    if( isolv==3 )
      ndof = 3*nel;   % OpenFOAM uses cell centered dofs.
    end
    data_i = [ data_i; [ isolv, ilev, nel, nvt, ndof, t_solv, it, err ] ];
    vel_profile = [ vel_profile; {p_r, p_z, w_vel, w_ref} ];

    if( iplot && ilev==nlev(end) )
      figure( 'name', solvers{isolv} )
      l_visualize( fea, solvers{isolv}, p_r, p_z, w_vel, w_ref, iimg )
    end

    fprintf( 1, '\n%s - Level %i - Error : %g\n\n', solvers{isolv}, ilev, err );

    delete(fullfile(ldir,lfile))
  end

  data = [ data; { data_i, vel_profile } ];
  if( isave )
    save 'ns8_benchmark_data' data
  end
end


% Data processing.
solvers = solvers(nsolv);
marker  = {'o','s','v'};
for i=1:length(solvers)

  data_i = data{i,1};

  fprintf('\n\n%s\n|------|---------|---------|------------|---------|-----|----------------|\n', solvers{i} )
  fprintf('| ilev |   nel   |   nvt   |    ndof    |  t_sol  |  it |      err       |\n' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|\n' )
  fprintf('| %4i | %7i | %7i | %10i | %7.1f | %3i | %14.11f |\n', data_i(:,2:end).' )
  fprintf('|------|---------|---------|------------|---------|-----|----------------|\n' )

  if( iplot )
    if( i==1 )
      figure
    end
    loglog( data_i(:,6), data_i(:,end), ['-',marker{i}], 'linewidth', 2 )
    hold on
    grid on
    xlabel( 'CPU time [s]' )
    ylabel( 'error' )
    title( 'Cost vs. accuracy' )
    if( i==length(data) )
      legend( solvers )
    end
    if( iimg )
      print('-r300','-dpng','ns8_benchmark_results');
    end
  end

end

if( ~nargout )
  clear data fea
end



%------------------------------------------------------------------------------%
function l_visualize( fea, solver, p_r, p_z, w_vel, w_ref, iimg )
% Postprocessing and visualization.

subplot(1,3,1)
postplot( fea, 'surfexpr', 'sqrt(u^2+w^2)', 'arrowexpr', {'u' 'w'} )
hold on
plot( p_r, p_z, 'k--' )
title( 'Velocity field' )

subplot(1,3,2)
postplot( fea, 'surfexpr', 'p', 'evaltype', 'exact', 'isoexpr', 'p' )
title( 'Pressure' )

subplot(1,3,3)
plot( p_r, w_ref, 'b-', 'linewidth', 2 )
hold on
plot( p_r, w_vel, 'r--', 'linewidth', 2 )
title('Velcity profile at outlet')
xlabel( 'Radius' )
legend( 'Reference', 'Computed', 'location', 'south' )
grid on

if( iimg )
  fname = ['ns8_',lower(solver),'_vis'];
  print('-r300','-dpng',fname);
end

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
