%FOAM2GRID Convert and import an OpenFOAM mesh.
%
%   [ GRID ] = GRID2FOAM( PATH, N_SDIM ) Imports an OpenFOAM mesh from
%   PATH (typically casedir/constant/polyMesh) and converts it to and
%   converts it to FEATool grid data.
%
%       boundary (only ysed for 2D meshes)
%       cellZones (for multiple subdomains, labeled ignored)
%       faces (list of vectors)
%       neighbour (unused)
%       owner (scalar list)
%       vertices (vector list)
%
%   (Note that meshes with mixed element types are not supported, and
%   boundary ordering and labels are ignored). N_SDIM optionally
%   specifies which mesh dimension to import (3 or 2 for 2D grids with
%   one mesh layer in the z-dimension).
%
%   See also GRID2FOAM

% Copyright 2013-2025 Precise Simulation, Ltd.
