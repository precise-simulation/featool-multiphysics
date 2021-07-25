function [ varargout ] = evalexpr0( varargin )
%EVALEXPR0 Evaluates an expression in a point on a group of cells.
%
%   [ VEVAL ] = EVALEXPR0( S_EXPR, XI, IND_S, IND_C, IND_E, PROB, AJAC ) Given a valid PROB
%   struct evaluates the expression S_EXPR on cells with index IND_C
%   at point XI in local coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 String expression to evaluate
%       xi          [n_sdim(+1),1]         Local coordinates of evaluation point
%       ind_s       scalar                 Subdomain to evaluate coefficients
%       ind_c       [n_c]                  Index vector to cells to evaluate
%       ind_e       scalar/array [n_c]     Index vector to edges to evaluate
%       prob        struct                 Problem definition struct
%                                          (uses prob.sol.u(:,end) for evaluations)
%       aJac        (n_c,:)                Optional array for Jac. inverse and determinant
%       gradrec     scalar                 Flag to use gradient recovery for derivatives
%                                          for first order Lagrange elements (default 0)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vEval       [n_c]                  Output vector of evaluated values
%
%   See also EVALDVAR, EVALSFUN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evalexpr0, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evalexpr0', varargin{:} );
if( ~nargout ), clear varargout; end
