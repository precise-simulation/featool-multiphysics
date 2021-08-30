function [ varargout ] = geom_loft_faces( varargin )
%GEOM_LOFT_FACES Loft faces to solid.
%
%   [ SOUT, NEW_TAG ] = GEOM_LOFT_FACES( SIN, TAGS, FACES ) Loft
%   selected FACES from geometry objects with TAGS to a solid geometry
%   object.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_loft_faces, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_loft_faces', varargin{:} );
if( ~nargout ), clear varargout; end
