function [ varargout ] = gobj_grid( varargin )
%GOBJ_GRID Create geometry object from grid.
%
%   [ GOBJ ] = GOBJ_GRID( GRID, TAG )
%   Creates a geometry object from a grid. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       grid        struct                    Grid struct
%       tag         string  {G1}              Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_grid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_grid', varargin{:} );
if( ~nargout ), clear varargout; end
