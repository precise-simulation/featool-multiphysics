function [ varargout ] = evalexprp( varargin )
%EVALEXPRP Evaluates an expression in the grid points.
%
%   [ VEVAL ] = EVALEXPRP( S_EXPR, PROB, SOLNUM, IND ) Evaluates the
%   expression S_EXPR in the points specified in PROB.GRID.P (in
%   global x/y/z-coordinates). PROB is a valid finite element problem
%   struct. Optional parameters are SOLNUM which indicates the chosen
%   solution, and IND for indexing which grid points are to be
%   evaluated for.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 String expression to evaluate
%       prob        struct                 Problem definition struct
%       solnum      scalar {n_sols}        Solution number/time to evaluate
%       ind         vector {1:n_p}         Index to grid points for which to evaluate
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vEval       [n_ind,1]              Output vector of evaluated values

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evalexprp, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evalexprp', varargin{:} );
if( ~nargout ), clear varargout; end
