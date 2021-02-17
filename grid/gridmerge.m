function [ varargout ] = gridmerge( varargin )
%GRIDMERGE Merge two grids.
%
%   [ GRID ] = GRIDMERGE( GRID1, IND1, GRID2, IND2, I_DELETE, I_CHECK ) Merges
%   the matching boundaries edges and faces specified by boundary numbers IND1
%   in GRID1 and IND2 in GRID2. The algorithm tries to match up the boundary
%   with the fewest boundary edges/faces and nodes to the closest corresponding
%   ones on GRID2. The optional argument I_DELETE selects which grid to delete
%   duplicate nodes from (I_DELETE is 2 per default). I_CHECK toggles error
%   checking on/off.
%
%   Can optionally be called as GRIDMERGE( GRID1, GRID2 ) where IND1 and IND2
%   are determined by any shared boundary numbers (also in this case boundaries
%   and subdomains are not renumbered).
%
%   Examples:
%
%      1) Create a 2D rivet by joining a rectangle and half circle.
%
%      grid1 = rectgrid( 4, 10, [0.3 0.7;0 1] );
%      grid2 = circgrid( 4, 3, 0.5, [0.5;1] );
%      grid2 = delcells( grid2, 'y<=1' );
%      grid  = gridmerge( grid1, 3, grid2, 4 );
%
%      2) Create a flow over cylinder benchmark grid by merging three grids.
%
%      grid1 = ringgrid( [0.05 0.06 0.08 0.11 0.15], 32, [], [], [0.2;0.2] );
%      grid2 = holegrid( 8, 1, [0 0.41;0 0.41], 0.15, [0.2;0.2] );
%      grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
%      grid3 = rectgrid( [0.41 0.5 0.7 1 1.4 1.8 2.2], 8, [0.41 2.2;0 0.41] );
%      grid  = gridmerge( grid3, 4, grid2, 6 );
%
%      3) Create a 3D grid with two brackets attached to an I-beam section.
%
%      grid01 = ringgrid( 1, 20, 0.03, 0.06, [0;0] );
%      indc01 = selcells( grid01, 'y<=sqrt(eps)' );
%      grid01 = delcells( grid01, indc01 );
%
%      grid02 = holegrid( 5, 1, .06*[-1 1;-1 1], .03, [0;0] );
%      indc02 = selcells( grid02, 'y>=-sqrt(eps)' );
%      grid02 = delcells( grid02, indc02 );
%      grid2d = gridmerge( grid01, [5 6], grid02, [7 8] );
%
%      grid1 = gridextrude( grid2d, 1, 0.02 );
%      grid1 = gridrotate( grid1, pi/2, 1 );
%      grid2 = grid1;
%      grid1.p(2,:) = grid1.p(2,:) + 0.03;
%      grid2.p(2,:) = grid2.p(2,:) - 0.01;
%
%      x_coord = [ -0.08 linspace(-0.06,0.06,6) 0.08 ];
%      y_coord = [ -0.2 -0.15 -0.1 -0.05 -0.03 -0.01 ...
%                   0.01  0.03  0.05  0.1  0.15  0.2 ];
%      grid3 = blockgrid( x_coord, y_coord, 1, ...
%                         [-0.08 0.08;-0.2 0.2;-0.08 -0.06] );
%      grid4 = blockgrid( 1, y_coord, 5, ...
%                         [-0.01 0.01;-0.2 0.2;-0.18 -0.08] );
%      grid5 = grid3;
%      grid5.p(3,:) = grid5.p(3,:) - 0.12;
%
%      grid = gridmerge( grid1, 8, grid3, 6 );
%      grid = gridmerge( grid2, 8, grid, 19 );
%      grid = gridmerge( grid4, 6, grid, 24, 1 );
%      grid = gridmerge( grid5, 6, grid, 33, 2 );
%
%   See also GRIDEXTRUDE, GRIDREVOLVE, GRIDROTATE, GRIDSCALE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridmerge, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridmerge', varargin{:} );
if( ~nargout ), clear varargout; end
