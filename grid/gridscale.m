function [ varargout ] = gridscale( varargin )
%GRIDSCALE Scale a grid.
%
%   [ GRID ] = GRIDSCALE( GRID, S ) Scaling a grid with the factors
%   specified in S. S can either be a real array with constant scaling
%   factors for each direction, for example S = [1 2 1] for a 2 times
%   scaling in the y-direction. Alternatively, S can be a cell array
%   with string expressions for the scaling factor.
%
%   Examples:
%
%      1) Scale a circle by a factor of 2 in the x-direction:
%
%      grid = spheregrid();
%      grid = gridscale( grid, [2 1 1] );
%
%      2) Scale the upper half of a rectangle to form an angle:
%
%      grid = rectgrid();
%      grid = gridscale( grid, {'1-(y>0.5).*(y-0.5)' 1} );
%
%   See also GRIDEXTRUDE, GRIDMERGE, GRIDREVOLVE, GRIDROTATE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridscale, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridscale', varargin{:} );
if( ~nargout ), clear varargout; end
