function [ varargout ] = geo2geom( varargin )
%GEO2GEOM Generate geometry objects from Gmsh geo file.
%
%   [ OBJS ] = GEO2GEOM( FILE ) Reads the Gmsh geo file, and
%   returns a cell array of geometry objects in OBJS.
%
%   [ P, E, F, V ] = GEO2GEOM( FILE ) Reads the Gmsh geo file,
%   extracts the points/vertices in P, edges E, faces F,
%   (and volumes in V).
%
%   See also GEOM2GEO

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geo2geom, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geo2geom', varargin{:} );
if( ~nargout ), clear varargout; end
