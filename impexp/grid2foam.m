function [ varargout ] = grid2foam( varargin )
%GRID2FOAM Generate/convert grid data to OpenFOAM format.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help grid2foam, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'grid2foam', varargin{:} );
if( ~nargout ), clear varargout; end
