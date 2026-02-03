function [ x, flag, t_solve ] = solvelin( A, b, type, x0, varargin )
%SOLVELIN Solve linear system Ax=b.
%
%   [ X, FLAG, T_SOLVE ] = SOLVELIN( A, B, TYPE, X0, VARARGIN ) Solves
%   the linear sparse system Ax = b with solver of TYPE (backslash,
%   mumps, gmres, bicgstab, or amg). X0 is an optional initial guess
%   for the iterative solver types (gmres/bicgstab/amg) and T_SOLVE is
%   the total time. FLAG returns 0 for success and ~0 otherwise.
%
%   See also SOLVESTAT, SOLVETIME

% Copyright 2013-2026 Precise Simulation, Ltd.

persistent try_mumps solvers

if( nargin<=1 && nargout==1 )  % Return available/valid solvers.
  if( isempty(solvers) )
    solvers = {'backslash', 'gmres', 'bicgstab'};

    ml2019b = verLessThan('matlab','9.8');
    use_mex = feacheckmex();

    if( use_mex && exist('dmumpsmex','file')==3  && ...
        ispc || (isunix && ~ismac) )
      solvers = [solvers, {'mumps'}];
    end

    if( use_mex && exist('amg','file')==3 && ml2019b )
      solvers = [solvers, {'amg'}];
    end

    solvers = sort(solvers);
  end

  if( nargin==1 )
    x = any(strcmpi(solvers,strtrim(A)));
  else
    x = solvers;
  end
  return
end

C_SOLVER_TYPES = {'backslash','mumps','gmres','bicgstab','amg'};

t0 = tic();
if( nargin<3 || isempty(type) )
  type = appconfig('linsolv');
  if( isempty(type) || (ischar(type) && ~any(strcmpi(type,C_SOLVER_TYPES))) )
    type = C_SOLVER_TYPES{1};
  end
end


flag = 0;
backslash = false;
if( ischar(type) )
  type = max( [1,find(strcmpi(type,C_SOLVER_TYPES))] );
end

if( ~(isscalar(type) && isnumeric(type) && type<=length(C_SOLVER_TYPES)) || ...
    (type==2 && isequal(try_mumps,false)) )
  type = 1;
end

try
  switch( type )

    case 2   % mumps
      if( ~(issparse(A) && isreal(A) && isreal(b)) )
        % warning( 'Complex valued systems not supported by mumps.' )
        backslash = true;
      else

        mumpsfunc = 0;
        if( length(b) < 500000 && isunix && ~ismac )
          mumpsfunc = 1;
        end

        % Initialization of a MATLAB Mumps structure.
        id = initmumps();
        id.SYM = 0;
        id = dmumps(id, [], mumpsfunc);
        cleanupObj1 = onCleanup(@()l_mumps_cleanup(id,mumpsfunc));
        assert(id.INFOG(1) == 0)
        id.JOB = 6;   % Set analysis + factorization + solve.
        id.ICNTL(1:4) = -1;   % suppress output.
        % id.ICNTL(1:4) = [6,0,6,2];   % standard output.

        % Mumps reordering:
        %    0 - Approximate Minimum Degree (AMD)
        %    2 - Approximate Minimuim Fill (AMF)
        %    3 - SCOTCH
        %    4 - PORD
        %    5 - METIS
        %    6 - Approximate Minimum Degree with automatic quasi row-detection (QAMD)
        %    7 - Automatic choice by MUMPS
        id.ICNTL(7) = 7;
        id.ICNTL(14) = 25;   % Percentage increase in estimated working space.
        id.RHS = b;   % Set RHS/load vector.

        % Call Mumps.
        id = dmumps(id, A, mumpsfunc);

        flag = id.INFOG(1);
        if( flag==0 )
          x = id.SOL;
        else
          warning( ['MUMPS linear solver failed with error INFO(1:2) = ', ...
                    num2str(flag),':',num2str(id.INFOG(2))] )
          x = zeros(size(b));
        end

        clear cleanupObj1
        t_solve = toc(t0);
      end

    case {3,4}   % gmres, bicgstab

      ILU_LEV = 4;
      RESTART = 100;
      TOL     = 1e-8;
      MAXIT   = 150;

      if( nargin<4 || isempty(x0) )
        x0 = zeros(size(b));
      end

      % Compute ILU preconditioner.
      if( ILU_LEV<=0 )
        ILU_SETUP.type    = 'crout';   % nofill, ilutp, crout
        ILU_SETUP.droptol = 1e-2;
        ILU_SETUP.milu    = 'row';     % row, col, off

        [L,U] = ilu(  A, ILU_SETUP );
      else
        [L,U] = iluk( A, ILU_LEV );
      end

      if( strcmp(lower(type),'gmres') )
        [x,flag,relres,iter,resvec] = gmres( A, b, RESTART, TOL, MAXIT, L, U, x0 );
      else
        [x,flag,relres,iter,resvec] = bicgstab( A, b, TOL, MAXIT, L, U, x0 );
      end
      if( flag~=0 || iter>=MAXIT )
        warning( [C_SOLVER_TYPES{type},' solver failed to converge.'] )
      end

      case 5   % amg(cl)

        if( nargin<4 || isempty(x0) )
          x0 = zeros(size(b));
        end

        cOptDef = {'tol', 1e-8;
                   'maxit', 150;
                   'reusemode', 1;
                   'preconditioner', 1;
                   'coarsening', 1;
                   'relaxation', 3;
                   's_relaxation', 3;
                   'solver', 4;
                   'block_size', 0;
                   'active_rows', 0;
                   'use_drs', 0;
                   'drs_eps_ps', 0.0200;
                   'drs_eps_dd', 0.2000;
                   'drs_row_weights', [];
                   'update_sprecond', 0;
                   'update_ptransfer', 0;
                   'cpr_blocksolver', 1;
                   'coarse_enough', -1;
                   'direct_coarse', 1;
                   'max_levels', -1;
                   'ncycle', 1;
                   'npre', 1;
                   'npost', 1;
                   'pre_cycles', 1;
                   'gmres_m', 30;
                   'lgmres_k', 3;
                   'lgmres_always_reset', 1;
                   'lgmres_store_av', 1;
                   'idrs_s', 4;
                   'idrs_omega', 0.7000;
                   'idrs_replacement', 0;
                   'bicgstabl_l', 2;
                   'bicgstabl_delta', 0;
                   'bicgstabl_convex', 1;
                   'aggr_eps_strong', 0.0800;
                   'aggr_over_interp', 1;
                   'rs_eps_strong', 0.2500;
                   'rs_trunc', 1;
                   'rs_eps_trunc', 0.2000;
                   'aggr_relax', 0.6667;
                   'write_params', 0;
                   'nthreads', 4;
                   'verbose', 0;
                   'ilut_p', 2;
                   's_ilut_p', 2;
                   'ilut_tau', 0.0100;
                   's_ilut_tau', 0.0100;
                   'iluk_k', 1;
                   's_iluk_k', 1;
                   'ilu_damping', 1;
                   's_ilu_damping', 1;
                   'jacobi_damping', 0.7200;
                   's_jacobi_damping', 0.7200;
                   'chebyshev_degree', 5;
                   's_chebyshev_degree', 5;
                   'chebyshev_lower', 0.0333;
                   's_chebyshev_lower', 0.0333;
                   'chebyshev_power_iters', 0;
                   's_chebyshev_power_iters', 0 };
        [got,val]  = parseopt(cOptDef,varargin{:});

        try
          id = 1;
          [ x, err, iter ] = amg( A.', b, val, val.tol, val.maxit, id, val.reusemode, x0 );
        catch
          x = x0;
          err = inf;
          iter = -1;
        end

        flag = err > val.tol;
        if( flag~=0 || iter>=val.maxit )
          warning( [C_SOLVER_TYPES{type},' solver failed to converge.'] )
        end

    otherwise
      backslash = true;
  end

catch
  warning( [C_SOLVER_TYPES{type},' linear solver failed. Reverting to built-in solver.'] )
  if( type==2 ), try_mumps = false; end
  backslash = true;
end

if( backslash )
  x = mldivide( A, b );
end

t_solve = toc(t0);


%------------------------------------------------------------------------------%
function l_mumps_cleanup(id, mumpsfunc)
% Release mumps memory.
id.JOB = -2;
id = dmumps(id, [], mumpsfunc);
