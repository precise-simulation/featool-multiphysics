function [ varargout ] = evalexprdof( varargin )
%EVALEXPRDOF Evaluates an expression in the degrees of freedom.
%
%   [ VEVAL, IND_DOF ] = EVALEXPRDOF( EXPR, DVAR_OR_SFUN, PROB, IND_B )
%   Given a valid PROB struct evaluates the expression(s) EXPR in the
%   degrees of freedom corresponding to dependent variable DVAR (or
%   shape function SFUN). IND_B indicates evaluation for selected boundaries
%   (no argument evaluates for all dofs).
%
%       Input       Value                  Description
%       -----------------------------------------------------------------------------------
%       expr        string (cell array)    String expression(s) to evaluate
%       dvar/sfun   scalar (or string)     Index to dvar or sfun string expression
%       prob        struct                 Problem definition struct
%       ind_b       integers               (Optional) evaluation boundaries
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vEval       [n_dof,n_expr]         Output vector/array of evaluated values
%       ind_dof     [n_dof,1]              Dof numbering for evaluated values
%
%   See also EVALEXPR0

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evalexprdof, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evalexprdof', varargin{:} );
if( ~nargout ), clear varargout; end
