function [ varargout ] = holegrid( varargin )
%HOLEGRID Generate 2d rectangular grid with a circular hole.
%
%   [ GRID ] = HOLEGRID( NS, NR, XP, R, XC, TH_OFFSET ) Generates a grid for
%   a rectangular domain with a circular hole. NS specifies the cell
%   resolution of the outer square (default 5), and NR the number of cells
%   in the radial direction of the outer layer (default 5). XP=[x1 y1 x2 y2]
%   specifies the coordinates of the outer lower left (x1,y1) and upper right
%   (x2,y2) corners.The optional arguments XC=[x_c0;y_c0] and R specify the
%   center coordinates and radius of the circle (default R=0.5 and XC=[0;0]).
%   Furthermore, TH_OFFSET specify a rotation of the whole grid.
%
%   Examples:
%
%      1) A 2x2 square with a circular hole with radius 0.5.
%
%      grid = holegrid( 5, 5, [-1 1;-1 1], 0.5, [0;0] );
%
%      2) Rotated grid with offset hole.
%
%      grid = holegrid( 5, 5, [-1 1;-1 1], 0.5, [0.25;0], pi/4 );
%
%   See also BLOCKGRID, CIRCGRID, CYLGRID, LINEGRID, RECTGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help holegrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'holegrid', varargin{:} );
if( ~nargout ), clear varargout; end
