function [ data, fea ] = run_featool_benchmarks( opt )
%RUN_FEATOOL_BENCHMARKS Run FEATool benchmark scripts.
%
%   [ DATA, FEA ] = RUN_FEATOOL_BENCHMARKS( OPT ) Run FEATool
%   benchmark scripts. Input is a struct OPT with the following fields
%
%       basename - benchmark case name
%       solvers  - cell array of solver name/descriptions
%       grids    - cell array of grid cell name/descriptions
%       cases    - cell array of benchmark case definitions (n_cases x 4)
%                      cases = { solver_type, grid_type, shape_functions, grid_levels; }
%       fcn_fea  - function handle for constructing fea data struct callable as
%                      fea = fcn_fea( grid_type, case, grid_level )
%       fcn_err  - function handle for computing errors, err = fcn_err( fea )
%       fcn_proc - (optional) function handle to post-process data
%
%   The function can also be called with a file string to process
%   previously saved benchmark data.
%
%   See also EX_NAVIERSTOKES1B, EX_NAVIERSTOKES2B, EX_NAVIERSTOKES3B, EX_NAVIERSTOKES8B

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ischar(opt) && exist(opt)==2 )
  load( opt );
  l_process_data( opt, data, fea );
  return
end

FEATOOL  = 1;
FENICS   = 2;
OPENFOAM = 3;
SU2      = 4;

% Required inputs.
basename = opt.basename;
solvers  = opt.solvers;
grids    = opt.grids;
cases    = opt.cases;
fcn_fea  = opt.fcn_fea;
fcn_err  = opt.fcn_err;

if( ~isfield(opt,'workdir') )   % Benchmark/output directory.
  opt.workdir = fullfile(pwd(),opt.basename);
end
rmkdir(opt.workdir);
logdir = fullfile(tempdir(),opt.basename);   % Logfile directory.
rmkdir(logdir);
logfile = [basename,'.log'];   % Logfile name.
logfullfile = fullfile(logdir,logfile);

diary(fullfile(opt.workdir,[opt.basename,'.txt']))
fprintf( 'START: %s %s\n\nSystem Info:\n\n', opt.workdir, datestr(now()) )
cleanupObj1 = onCleanup(@()cleanup(logdir));
cleanupObj2 = onCleanup(@()diary('off'));

try
  [~,sys,os,cpu,ml] = sysinfo();
  fprintf( '  %s\n', [sys,' ',cpu], os, ml );
catch,end

% Run benchmarks.
data = {};
for i_case=1:size(cases,1)
  case_i = cases(i_case,:);
  i_solver = case_i{1};
  i_grid   = case_i{2};
  sfun_i   = case_i{3};
  for j=1:length(sfun_i)
    ipos_ = max([1, find( sfun_i{j} == '_', 1, 'last' ) + 1]);
    sfun_i{j} = strrep( sfun_i{j}(ipos_:end), 'disc0', 'P0' );
    sfun_i{j} = strrep( sfun_i{j}, 'disc1', 'P-1' );
  end
  s_case = [solvers{i_solver},'-',grids{i_grid},'-',strcat(sfun_i{:})];
  solver_args = {};
  if( length(case_i)>=5 )
    solver_args = case_i{5};
  end


  data_i = [];
  for i_lev=case_i{4}

    fea = fcn_fea( i_grid, case_i{3}, i_lev );

    fprintf( '\n%s - Level %i\n\n', s_case, i_lev );

    % Parse and solve problem.
    fid = [];
    try
      switch( i_solver )

        case {FEATOOL}
          fid = fopen( logfullfile, 'w+' );
          fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 1, 'tolchg', 1e-4, 'relchg', 0, 'maxnit', 50, solver_args{:} );

        case {FENICS}     % Use external FEniCS solver.
          fea = fenics( fea, 'fdir', logdir, 'fname', opt.basename, 'clean', false, 'maxnit', 50, solver_args{:} );

        case {OPENFOAM}   % Use external OpenFOAM CFD solver.
          fea.sol.u = openfoam( fea, 'casedir', logdir, 'logfname', logfile, 'clean', false, ...
                                'interp', 2*(~strcmp(fea.sfun{1},'sf_disc0')), ...
                                'upwind', 0, 'tolres', 1e-6, 'endTime', 2000, 'writePrecision', 16, solver_args{:} );

        case {SU2}        % Use external SU2 CFD solver.
          fea.sol.u = su2( fea, 'workdir', logdir, 'logfname', logfile, 'clean', false, solver_args{:} );
      end

      [t_sol,it] = l_parse_logfile( i_solver, logfullfile, fid );

      % Error checking.
      err = fcn_err( fea )

    catch me
      warning(me.message)
      t_sol = nan;
      it    = nan;
      err   = nan;
    end

    nel = size(fea.grid.c,2);
    nvt = size(fea.grid.p,2);
    if( i_solver==OPENFOAM && strcmp(sfun_i{1},'P0') )
      ndof = length(fea.dvar)*nel;   % If OpenFOAM uses cell centered dofs.
    else
      ndof = sum(fea.eqn.ndof);
    end
    data_i = [ data_i; [ i_lev, nel, nvt, ndof, t_sol, it, err ] ];

    cleanup( logfile )
  end

  data = [ data; { s_case, data_i } ];

  savefile = fullfile(opt.workdir,opt.basename);
  save( savefile, 'fea', 'data', 'opt' );
end
fprintf( '\nEND: %s %s\n\n', opt.workdir, datestr(now()) )


% Data processing.
l_process_data( opt, data, fea );


%------------------------------------------------------------------------------%
function l_process_data( opt, data, fea )

if( isfield(opt,'fcn_proc') && isa(opt.fcn_proc,'function_handle') )
  opt.fcn_proc( opt, data, fea );
  return
end

COLORS  = {'r','g','b','c','m','y',[.7 .7 1],[.7 1 .7], ...
           [1 .7 .7],[1 0.5 0],[0.5 1 0],[1 0 0.5],'k'};
MARKERS = {'o','s','v','^','+','x','d','<','>','p','h','*'};

h_fig = figure(); hold on
for i=1:size(data,1)

  s_case = data{i,1};
  data_i = data{i,2};

  fprintf( '\n\n%s\n|------|---------|---------|------------|---------|-----|----------------|\n', s_case )
  fprintf( '| ilev |   nel   |   nvt   |    ndof    |  t_sol  |  it |      err       |\n' )
  fprintf( '|------|---------|---------|------------|---------|-----|----------------|\n' )
  fprintf( '| %4i | %7i | %7i | %10i | %7.1f | %3i | %14.11f |\n', data_i.' )
  fprintf( '|------|---------|---------|------------|---------|-----|----------------|\n' )

  h(i) = plot( data_i(:,5), abs(data_i(:,end)), ...
               ['-',MARKERS{mod(i-1,length(MARKERS))+1}], ...
               'color', COLORS{mod(i-1,length(COLORS))+1}, ...
               'linewidth', 2 );
end
set( gca, 'XScale', 'log', 'YScale', 'log' )

grid on
xlabel( 'CPU Time [s]' )
ylabel( 'Error' )
title(  'Cost vs. Accuracy' )

legend( data(:,1), 'Location', 'best' )

imgfile = fullfile(opt.workdir,opt.basename);
print( '-r300', '-djpeg', imgfile )
print( '-r300', '-dpng',  imgfile )

%------------------------------------------------------------------------------%
function [ t_solv, it ] = l_parse_logfile( i_solver, logfile, fid )
% Parse log file to determine solver time and number of iterations.

switch( i_solver )

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

    if( any(fopen('all')==fid) )
      fclose( fid );
    end

  case 2   % FEniCS
    fid = fopen( logfile, 'r' );
    cleanupObj = onCleanup( @()fclose(fid) );
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
    fid = fopen( logfile, 'r' );
    cleanupObj = onCleanup( @()fclose(fid) );
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

  case 4   % SU2
    fid = fopen( logfile, 'r' );
    cleanupObj = onCleanup( @()fclose(fid) );
    it = nan;
    while 1
      sline = fgetl( fid );
      if ~ischar(sline), break, end
      if isempty(sline), continue, end
      if( any(findstr(sline,'SU2 solution time:')) )
        t_solv = str2num(sline(19:end));
      else
        c = regexp( strtrim(sline), '\|', 'split' );
        c = c(~cellfun(@isempty,c));
        if( ~isempty(str2num(c{1})) )
          it = str2num(c{1});
        end
      end
    end

end
