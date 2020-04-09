function [ x, flag, t_solve ] = solvelin( A, b, type, x0 )
%SOLVELIN Solve linear system Ax=b.
%
%   [ X, FLAG, T_SOLVE ] = SOLVELIN( A, B, TYPE, X0 ) Solves the
%   linear sprase system Ax = b with solver TYPE (backslash, mumps,
%   gmres, or bicgstab). X0 is an optional initial guess for the
%   iterative solver types (gmres/bicgstab) and T_SOLVE is the total
%   time. The FLAG returns 0 for success and >0 otherwise.
%
%   See also SOLVESTAT, SOLVETIME

% Copyright 2013-2020 Precise Simulation, Ltd.


t0 = tic();
if( nargin<3 )
  type = 'backslash';
end


flag = 0;
backslash = false;
if( ischar(type) )
  c_solver_types = {'backslash','mumps','gmres','bicgstab'};
  type = max( [1,find(strcmpi(type,c_solver_types))] );
end
try
  switch( lower(type) )

    case 2   % mumps
      if( ~isreal(A) || ~isreal(b) )
        % warning( 'Complex valued systems not supported by mumps.' )
        backslash = true;
      else

        % Initialization of a MATLAB Mumps structure.
        id = initmumps();
        id.SYM = 0;
        id = dmumps(id);
        id.JOB = 6;   % Set analysis + factorization + solve.
        id.ICNTL(1:4) = -1;   %  Suppress output.

        % Mumps reordering:
        %    0 - Approximate Minimum Degree is used
        %    3 - SCOTCH (not available)
        %    4 - PORD (if available)
        %    5 - METIS (if available)
        %    7 - Automatic choice by MUMPS
        id.ICNTL(7) = 7;
        id.RHS = b;   % Set RHS/load vector.

        % Call Mumps.
        id = dmumps(id,A);
        x = id.SOL;

        % Release memory.
        id.JOB = -2;
        id = dmumps(id);
        t_solve = toc(t0);
      end

    case {3,4}   % gmres,bicgstab

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
      if( flag>0 )
        warning( [type,' solver failed to converge.'] )
      end

    otherwise
      backslash = true;
  end

catch
  warning( ['Linear solver ',type,' failed. Reverting to built-in solver.'] )
  backslash = true;
end

if( backslash )
  x = mldivide( A, b );
end

t_solve = toc(t0);
