function [ varargout ] = cylgrid( varargin )
%CYLGRID Generate 3d cylindrical hexahedral grid.
%
%   [ GRID ] = CYLGRID( NS, NR, NZ, R, LZ, XP, AX ) Generates a hexahedral
%   grid for a cylindrical domain with NS+NR cells in the radial direction.
%   NS specifies the cell resolution of the inner square (default 4), and NR
%   ns number of cells in the radial direction of the outer layer (default 3).
%   NZ specity the numver of cells in the lengthwise direction and LZ the
%   corresponding length. The optional arguments R, XP, and AX specify the
%   radius, center coordinates of the cylinder (default R = 1 and XP = [0;0;0]),
%   and the axis of alignment (default 1 equals to the x-axis).
%
%   Examples:
%
%      1) A cylindrical grid with radius 1 and length 1:
%
%      grid = cylgrid();
%
%      2) Cylinder with radius 0.5 length 2 centered and extending
%         in the negative y-direction from the point [1 1 2].
%
%      grid = cylgrid( 3, 4, 10, 0.5, 2, [1;1;2], -2 );
%
%   See also BLOCKGRID, CIRCGRID, HOLEGRID, LINEGRID, RECTGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help cylgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'cylgrid', varargin{:} );
if( ~nargout ), clear varargout; end
