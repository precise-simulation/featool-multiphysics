function [ varargout ] = evalexpr( varargin )
%EVALEXPR Evaluates an expression in specified points.
%
%   [ VEVAL ] = EVALEXPR( S_EXPR, XP, PROB, SOLNUM, EVAL_TYPE, TOL )
%   Evaluates the expression S_EXPR in the points specified in XP (in
%   global x/y/z-coordinates). PROB is a valid finite element problem
%   struct. The EVAL_TYPE flag (default 1) toggles vectorized
%   evaluation for dof variables with P0-P2 fem shape functions, or
%   allows faster evalution via intermediate nodal linear
%   interpolation. TOL is a tolerance used by ismembertol/deduplicate
%   to determined which if any evaluation points in XP are aligned
%   with grid points.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 String expression to evaluate
%       xp          [n_sdim,n_xp]          Coordinates of evaluation points
%       prob        struct                 Problem definition struct
%       solnum      scalar  {n_sols}       Solution number/time to evaluate
%       eval_type   scalar  {1}            Evaluation type
%                                             0 - Standard/full evaluation
%                                             1 - Vectorized evaluation for dep. variables
%                                             2 - Linear interpolation via nodes (fast)
%                                             3 - Like 2 and allow invalid (nan) values
%       tol         scalar  {1e-5}         Tolerance for evalation/grid points
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vEval       [n_xp,1]               Output vector of evaluated values

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evalexpr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evalexpr', varargin{:} );
if( ~nargout ), clear varargout; end
