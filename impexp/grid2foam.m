%GRID2FOAM Convert and export grid data to OpenFOAM format.
%
%   [ POLYMESH ] = GRID2FOAM( GRID, PATH, FID_LOG ) Converts and
%   exports FEATool grid data to the OpenFOAM ASCII mesh format.
%
%   Input data is a valid FEA or GRID struct, and a PATH to where the
%   OpenFOAM mesh files should be written (defaults to the current
%   directory). Exports the following files:
%
%       boundary (dictionary of boundary types, labeled boundary 1-n)
%       cellZones (for multiple subdomains, labeled zone1-n)
%       faces (list of vectors)
%       neighbour (scalar list)
%       owner (scalar list)
%       vertices (vector list)
%
%   FID_LOG optionally specifies a message log file handle (negative
%   for gui output or empty for no output).
%
%   Optionally returns grid dictionary definitions as a struct POLYMESH.
%
%   See also FOAM2GRID, OPENFOAM_DICT_POLYMESH, OPENFOAM_WRITE

% Copyright 2013-2025 Precise Simulation, Ltd.
