function [ varargout ] = solveeig( varargin )
%SOLVEEIG Solve eigenvalue problem.
%
%   [ U, L ] = SOLVEEIG( PROB, VARARGIN ) Solves the eigenvalue
%   problem described in the PROB finite element struct with
%   homogenous boundary conditions, and returns the eigenvectors U and
%   eigenvalues L. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       neigs       scalar/{6}                   Number of eigenvalues
%       sigma       scalar/string/{sm}           Range or type of eigenvalues
%                                                  see the eigs function for options
%       imass       scalar/{4}                   Mass matrix lumping
%                                                  1 - Full mass matrix
%                                                  2 - row sum lumping
%                                                  3 - diagonal lumping
%                                                  4 - HRZ diagonal lumping
%       init        u0|{expr}/{0}                Initial value solution or expression
%                                                for linearization point
%       solcomp     {all dvars}                  Dep. variables/subdomains to solve for
%       icub        scalar/{auto}                Cubature rule/order used in assembly
%                                                Default 1+max(shape function order)
%       isymm       scalar/{0}                   Symmetrize BCs if applicable
%       waitbar     scalar/{0}                   Show waitbar
%       fid         scalar/{1}                   File identifier for output ([]=no output)
%
%   See also EIGS, SOLVESTAT, SOLVETIME

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help solveeig, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'solveeig', varargin{:} );
if( ~nargout ), clear varargout; end
