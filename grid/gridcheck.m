function [ varargout ] = gridcheck( varargin )
%GRIDCHECK Check grid for errors.
%
%   [ N_ERR, CIND ] = GRIDCHECK( G, ISWARN, I_CHECK ) Checks a grid G for errors.
%   N_ERR returns the number of errors found. CIND is a cell array giving the
%   indices to incorrectly oriented cells, collapsed cells (with zero volume/area),
%   and non-convex cells. ISWARN is a flag to enable/disable warnings. I_CHECK is
%   a flag to choose check/test type (1 = test Jacobian determinant (default),
%   2 = test if polygon points are in clockwise order, 3 = split hexahedra into
%   six tetrahedra).
%
%   See also REORIENT_CELLS

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridcheck, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridcheck', varargin{:} );
if( ~nargout ), clear varargout; end
