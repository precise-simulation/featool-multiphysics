function [ varargout ] = copy_geometry_object( varargin )
%COPY_GEOMETRY_OBJECT Copy geometry object(s) with offset.
%
%   [ GEOM ] = COPY_GEOMETRY_OBJECT( TAGS, GEOM, OFFSET ) Copies the
%   geometry objects in GEOM with corresponding TAGS. An translational
%   OFFSET (default 0) is optional.

% Copyright 2013-2019 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help copy_geometry_object, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'copy_geometry_object', varargin{:} );
if( ~nargout ), clear varargout; end
