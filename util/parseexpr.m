function [ varargout ] = parseexpr( varargin )
%PARSEEXPR Parses a string expression.
%
%   [ S_EXPR, CSUBEXPR ] = PARSEEXPR( S_EXPR, S_SUBEXPR, C_ALLOWEXPR )
%   Parses the string expression in S_EXPR and returns a cell array
%   CSUBEXPR with variable expressions to process. Also returns a modified
%   expression S_EXPR which has the sub variable names substituted with
%   'S_SUBEXPR(:,I)' where I is the index in CSUBEXPR to the corresponding
%   variable. This enables a call to EVAL(S_EXPR) after the expressions in
%   CSUBEXPR has been computed and copied to the S_SUBEXPR array. The
%   optional array C_ALLOWEXPR allows strings as subexpressions (which
%   function names would otherwise be removed).
%
%   See also SYMVARBE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help parseexpr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'parseexpr', varargin{:} );
if( ~nargout ), clear varargout; end
