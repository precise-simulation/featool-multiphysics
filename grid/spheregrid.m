function [ varargout ] = spheregrid( varargin )
%SPHEREGRID Generate 3d hexahedral grid for a sphere.
%
%   [ GRID ] = SPHEREGRID( NS, NR, R, XP ) Generates a hexahedral grid for a
%   spherical domain with NS+NR cells in the radial direction. NS specifies
%   the cell resolution of the inner cube (default 4), and NR the number of
%   cells in the radial direction of the outer layer (default 3). The optional
%   arguments R and XP = [x0;y0;z0] specify the radius and center coordinates
%   of the sphere (default R = 1 and XP = [0;0;0]).
%
%   Examples:
%
%      1) Generate a grid for a sphere with radius 1 centered at [ 0, 0 ]:
%
%      grid = spheregrid();
%
%      2) Generate a grid for a sphere with radius 0.5 centered at [ 1, 1, 1 ]:
%
%      grid = spheregrid( 16, 12, 0.5, [1;1;1] );
%
%   See also BLOCKGRID, CIRCGRID, CYLGRID, HOLEGRID, LINEGRID, RECTGRID, RINGGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help spheregrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'spheregrid', varargin{:} );
if( ~nargout ), clear varargout; end
