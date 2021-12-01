function [ varargout ] = gridrevolve( varargin )
%GRIDREVOLVE Revolve and extrude a 2D grid.
%
%   [ GRID ] = GRIDREVOLVE( GRID, N_CR, LZ, N_REV, R_OFF, TH_XZ, I_JOIN )
%   Revolves and extrudes a 2D GRID to generate a 3D grid. N_CR is the
%   number of cell rows in the revolved direction, and LZ the
%   length/offset height in the z-direction. N_REV specifies the
%   number of revolutions, and R_OFF the radial offset of the original
%   input grid. TH_XZ optionally prescribes an initial angle of the
%   input grid plane (default pi/2). (I_JOIN is optionally a flag to
%   enforce joining the two ends, default I_JOIN = LZ==0 & N_REV==1).
%
%   Examples:
%
%      1) Create a torus grid:
%
%      grid = gridrevolve( circgrid(), 20, 0, 1, 2 );
%
%      2) Create a spiral square grid:
%
%      grid = gridrevolve( rectgrid(3), 100, 10, 3, 2 );
%
%      3) Create a cored apple grid:
%
%      grid = circgrid();
%      grid = delcells( grid, 'x>=0' );
%      grid = gridrevolve( grid, 20, 0, 1, 1 );
%
%   See also GRIDEXTRUDE, GRIDMERGE, GRIDROTATE, GRIDSCALE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridrevolve, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridrevolve', varargin{:} );
if( ~nargout ), clear varargout; end
