function [ varargout ] = evalinit( varargin )
%EVALINIT Calculates initial solution vector from given expressions.
%
%   [ U0 ] = EVALINIT( PROB, C_EXPR ) Given a valid PROB struct
%   evaluates the inital value expressions in C_EXPR (or PHYS struct
%   if C_EXPR is not given), and returns a solution vector U0.
%
%   See also EVALEXPR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evalinit, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evalinit', varargin{:} );
if( ~nargout ), clear varargout; end
