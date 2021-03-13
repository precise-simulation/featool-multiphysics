function [ varargout ] = ringgrid( varargin )
%RINGGRID Generate 2d grid of a ring.
%
%   [ GRID ] = RINGGRID( NR, NTH, RI, RO, XP, TH_START, TH_END )
%   Generates a quadrilateral grid for a ring shaped domain with NR
%   cells in the radial and NTH cells in the tangential directions,
%   respectively. RI and RO give the inner and outer radius, while
%   the optional arguments XP, TH_START, and TH_END specify the center
%   coordinates start, and end angles. Furthermore, the position
%   of the radial layers can be manually specified through a vector
%   of coordinates in NR (which overrides RI and RO).
%
%   Examples:
%
%      1) An annulus 64 cell grid for a circle with radius 1 and
%         outer radius 3 centered at [ 0.5, 0.5 ]:
%
%      grid = ringgrid( 10, 36, 1, 3, [0.5;0.5] );
%
%      2) An annulus with prescribed radial layers.
%
%      grid = ringgrid( [0.05 0.06 0.08 0.11 0.15], 32, [], [], [0;0], -3/4*pi );
%
%      3) Quarter ring.
%
%      grid = ringgrid( [0.05 0.06 0.08 0.11 0.15], 32, [], [], [0;0], 0, pi/2 );
%
%   See also BLOCKGRID, CIRCGRID, CYLGRID, HOLEGRID, LINEGRID, RECTGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help ringgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'ringgrid', varargin{:} );
if( ~nargout ), clear varargout; end
