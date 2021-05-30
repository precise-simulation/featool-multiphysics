function [ varargout ] = solvestat( varargin )
%SOLVESTAT Solve stationary problem.
%
%   [ U ] = SOLVESTAT( PROB, VARARGIN ) Solves the stationary problem described in
%   the PROB finite element struct. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       linsolv     string/{mumps}               Sparse linear solver
%                                                    backslash, mumps, gmres, bicgstab, amg
%       nsolve      scalar/{1}                   Non-linear solution method
%                                                  1 - Fixed point (Picard) iteration
%                                                  2 - Newton method
%                                                 <0 - |nsolve|-Picard steps before Newton
%       jac         scalar/struct {0}            Newton Jacobian computation type
%                                                 struct - direct assembly using symbolic
%                                                   jac.form and jac.coef struct contents
%                                                  0 - using numerical differentiation
%       minnit      scalar/{0}                   Minimum number of non-linear iterations
%       maxnit      scalar/{25}                  Maximum number of non-linear iterations
%       nlrlx       scalar/string/{1.0}          Relaxation between non-linear iterations,
%                                                scalar or string expr. ex 0.5*(1+(it>5))
%       nlinasm     logical/{[0 0 0]}            Force non-linear reassembly of A, f, bdrn
%       toldef      scalar/{1e-6}                Stopping criteria for solution defects
%       tolchg      scalar/{1e-6}                Stopping criteria for solution changes
%       reldef      logical/{0}                  Check relative defect changes
%       relchg      logical/{1}                  Check relative solution changes
%       init        u0|{expr}/{0}                Initial value solution or expression
%       solcomp     {all dvars}                  Dep. variables/subdomains to solve for
%       icub        scalar/{auto}                Cubature rule/order used in assembly
%                                                Default 1+max(shape function order)
%       isymm       scalar/{0}                   Symmetrize BCs if applicable
%       waitbar     scalar/{0}                   Show waitbar
%       fid         scalar/{1}                   File identifier for output ([]=no output)
%
%   See also SOLVETIME, SOLVEEIG

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help solvestat, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'solvestat', varargin{:} );
if( ~nargout ), clear varargout; end
