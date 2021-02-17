function [ varargout ] = linegrid( varargin )
%LINEGRID Generate 1D line grid.
%
%   [ GRID ] = LINEGRID( P, XMIN, XMAX ) Generates a regular 1D line grid with
%   coordinates in P. If P is an integer scalar {P=10 per default} P
%   equally spaced cells will be generated on a line spanning XMIN and XMAX.
%   If the optional arguments (only used if P is a scalar) XMIN and XMAX are
%   omitted they default to 0 and 1 respectively.
%
%   Examples:
%
%      1) 10 grid cells on a line spanned by [ 0..1 ]:
%
%      grid = linegrid();
%
%      2) 25 grid cells on a line spanned by [ -0.5..0.5 ]:
%
%      grid = linegrid( 25, -0.5, 0.5 );
%
%      3) 10 grid cells on a log spaced line [ -2..2 ]:
%
%      p = logspace(-2,2,11);
%      grid = linegrid( p );
%
%   See also BLOCKGRID, CIRCGRID, CYLGRID, HOLEGRID, RECTGRID, RINGGRID, SPHEREGRID

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help linegrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'linegrid', varargin{:} );
if( ~nargout ), clear varargout; end
