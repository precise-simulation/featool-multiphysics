function [ varargout ] = circgrid( varargin )
%CIRCGRID Generate 2d quadrilateral grid for a circle.
%
%   [ GRID ] = CIRCGRID( NS, NR, R, XP, TH_OFFSET ) Generates a quadrilateral
%   grid for a circular domain with NS+NR cells in the radial direction. NS
%   specifies the cell resolution of the inner square (default 4), and NR the
%   number of cells in the radial direction of the outer layer (default 3).
%   The optional arguments R and XP = [x0;y0] specify the radius and center
%   coordinates of the circle (default R = 1 and XP = [0;0]). Furthermore,
%   TH_OFFSET specifies a rotation of the whole grid in radians.
%
%   Examples:
%
%      1) A 64 cell grid for a circle with radius 1 centered at [ 0, 0 ]:
%
%      grid = circgrid();
%
%      2) A 1024 cell grid for a circle with radius 0.5 centered at [ 1, 1 ]:
%
%      grid = circgrid( 16, 12, 0.5, [1;1], 0 );
%
%   See also CYLGRID, BLOCKGRID, HOLEGRID, LINEGRID, RECTGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help circgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'circgrid', varargin{:} );
if( ~nargout ), clear varargout; end
