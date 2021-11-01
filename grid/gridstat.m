function [ varargout ] = gridstat( varargin )
%GRIDSTAT Grid statistics.
%
%   [ STAT, S ] = GRIDSTAT( GRID ) Computes and returns grid
%   statistics, such as the grid cell area and volume, and the quality
%   measure; radius ratio for simplices and min/max diagonal ratio for
%   quadrilaterals and hexahedrons.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridstat, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridstat', varargin{:} );
if( ~nargout ), clear varargout; end
