function [ varargout ] = blockgrid( varargin )
%BLOCKGRID Generate 3D regular tensor product hexahedral grid.
%
%   [ GRID ] = BLOCKGRID( N_CX, N_CY, N_CZ, XP ) Generates a tensor
%   product block grid with N_CX, N_CY, NCZ hexahedral cells in the x,
%   y, and z-directions. The optional argument XP=[x1 x2; y1 y2; z1
%   z2] specifies the coordinates of the lower left (x1,y1,z1) and
%   upper right (x2,y2,z2) corners. N_CX, N_CY, and N_CZ can also be
%   vectors of x, y and z-coordinates. A call without input arguments
%   generates a default 10x10x10 unit cube.
%
%   Examples:
%
%      1) A 10x10x10 grid on the domain spanned by [ 0..1; 0..1; 0..1 ]:
%
%      grid = blockgrid();
%
%      2) A 20x40x80 grid on the domain spanned by [ -0.5..0.5; 0..2; 0..1 ]:
%
%      grid = blockgrid( 20, 40, 80, [-0.5 0.5; 1 2; 0 1] );
%
%      3) A grid spanning [0, 0.1, 0.4, 1], [0:0.25:1.5], [0, 1]:
%
%      grid = blockgrid( [0, 0.1, 0.4, 1], [0:0.25:1.5], [0, 1] );
%
%   See also CIRCGRID, CYLGRID, HOLEGRID, LINEGRID, RECTGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help blockgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'blockgrid', varargin{:} );
if( ~nargout ), clear varargout; end
