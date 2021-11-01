function [ varargout ] = evaldvar( varargin )
%EVALDVAR Evaluates a dependent variable in a point on a group of cells.
%
%   [ VDVARVAL ] = EVALDVAR( SFUN, I_EVAL, N_SDIM, N_VERT, XI, AINVJAC, DOFMAP, U )
%   Evaluates a dependent variable in point XI on the cells given in dofmap.
%   SFUN specifies which shape function to use with the corresponding
%   solution vector U.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       sfun        string                 Shape function string
%       i_eval      scalar                 Evaluation flag
%       n_sdim      scalar                 Number of space dimensions
%       n_vert      scalar                 Number of vertices per cell
%       xi          [n_sdim(+1)]           Local coordinates of evaluation point
%       aInvJac     [n,n_sdim(+1)*n_sdim]  Inverse of transformation Jacobian
%       dofmap      [n_ldof,n_cells]       Dof connectivity map (for each cell)
%       u           [n_u]                  Solution vector
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vDVarVal    [n_cells]              Evaluated values
%
%   See also EVALEXPR, EVALSFUN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help evaldvar, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'evaldvar', varargin{:} );
if( ~nargout ), clear varargout; end
