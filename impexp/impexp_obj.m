%IMPEXP_OBJ Import/export geometry in OBJ format.
%
%   [ GEOM ] = IMPEXP_OBJ( FILE_NAME, MODE, GEOM, TEST2D ) Import and
%   export of geometries in Wavefront OBJ format. Triangular facets
%   with the fields v for vertices, f for faces are supported.
%
%   FILE_NAME is a string which specifies the file to process.
%
%   MODE is a string indicating either IMPORT (default) or EXPORT.
%
%   During export, a geometry GEOM must be given. If GEOM is provided
%   during import the output will include the original geometry objects.

% Copyright 2013-2022 Precise Simulation, Ltd.
