%SOLVELIN Solve linear system Ax=b.
%
%   [ X, FLAG, T_SOLVE ] = SOLVELIN( A, B, TYPE, X0, VARARGIN ) Solves
%   the linear sparse system Ax = b with solver of TYPE (backslash,
%   mumps, gmres, bicgstab, or amg). X0 is an optional initial guess
%   for the iterative solver types (gmres/bicgstab/amg) and T_SOLVE is
%   the total time. FLAG returns 0 for success and ~0 otherwise.
%
%   See also SOLVESTAT, SOLVETIME

% Copyright 2013-2025 Precise Simulation, Ltd.

































