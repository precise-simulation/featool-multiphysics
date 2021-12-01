function [ varargout ] = solvetime( varargin )
%SOLVETIME Solve time-dependent problem.
%
%   [ U, TLIST ] = SOLVETIME( PROB, VARARGIN ) Solves the time-dependent problem described
%   in the PROB finite element struct. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       linsolv     string/{mumps}               Sparse linear solver
%                                                    backslash, mumps, gmres, bicgstab, amg
%       tstep       scalar/{0.1}                 Time step size (average for FS-scheme)
%       tmax        scalar/{1}                   Maximum simulation time
%       tstop       scalar/{1e-6}                Stopping criterita
%                                                (for solution changes in time)
%       ischeme     scalar/{2}                   Time stepping scheme
%                                                  1 - 1st order Backward Euler
%                                                  2 - 2nd order Crank-Nicolson
%                                                  3 - 2nd order Fractional-step-theta
%       imass       scalar/{4}                   Mass matrix lumping
%                                                  1 - Full mass matrix
%                                                  2 - row sum lumping
%                                                  3 - diagonal lumping
%                                                  4 - HRZ diagonal lumping
%       icub        scalar/{auto}                Cubature rule/order used in assembly
%                                                Default 1+max(shape function order)
%       minnit      scalar/{0}                   Minimum number of non-linear iterations
%       maxnit      scalar/{20}                  Maximum number of non-linear iterations
%       nstbwe      scalar/{0}                   Number of initial forced BE steps
%       nlrlx       scalar/string/{1.0}          Relaxation for non-linear iter., scalar
%                                                or string expression, ex (1+(t_sim>5))/2
%       nlinasm     logical/{[0 0 0 0]}          Force reassembly of M, A, f, bdrn
%       toldef      scalar/{1e-6}                Stopping criteria for solution defects
%       tolchg      scalar/{1e-6}                Stopping criteria for solution changes
%       reldef      logical/{0}                  Check relative defect changes
%       relchg      logical/{0}                  Check relative solution changes
%       init        u0|{expr}/{0}                Initial value solution or expression
%       solcomp     {all dvars}                  Dep. variables/subdomains to solve for
%       waitbar     scalar/{0}                   Show waitbar
%       fid         scalar/{1}                   File identifier for output ([]=no output)
%
%   See also SOLVESTAT, SOLVEEIG

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help solvetime, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'solvetime', varargin{:} );
if( ~nargout ), clear varargout; end
