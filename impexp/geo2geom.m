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

% Copyright 2013-2025 Precise Simulation, Ltd.
