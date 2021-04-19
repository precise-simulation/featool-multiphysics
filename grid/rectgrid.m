function [ varargout ] = rectgrid( varargin )
%RECTGRID Generate 2D rectangular tensor product quadrilateral grid.
%
%   [ GRID ] = RECTGRID( N_CX, N_CY, XP ) Generates a tensor product
%   quadrilateral grid with N_CX, N_CY cells in the x and
%   y-directions, respectively. The optional argument XP=[x1 x2; y1
%   y2] specifies the coordinates of the lower left (x1,y1) and upper
%   right (x2,y2) corners. N_CX and N_CY can also be vectors of x and
%   y-coordinates. A call without input arguments generates a default
%   10x10 unit square.
%
%   Examples:
%
%      1) A 10x10 grid on the domain spanned by [ 0..1; 0..1 ]:
%
%      grid = rectgrid();
%
%      2) A 20x40 grid on the domain spanned by [ -0.5..0.5; 1..2 ]:
%
%      grid = rectgrid( 20, 40, [-0.5 0.5; 1 2] );
%
%      3) A grid spanning [0, 0.1, 0.4, 1], [0:0.25:1.5]:
%
%      grid = rectgrid( [0, 0.1, 0.4, 1], [0:0.25:1.5] );
%
%   See also BLOCKGRID, CIRCGRID, CYLGRID, HOLEGRID, LINEGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help rectgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'rectgrid', varargin{:} );
if( ~nargout ), clear varargout; end
