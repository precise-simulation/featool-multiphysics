%GEOM2POLY Generate poly data for use with Triangle.
%
%   [ POLY ] = GEOM2POLY( GEOM, HMAX, HMAX_E, FILE_NAME, FID_LOG ) Generate poly data
%   for use with the 2D unstructured grid generator Triangle. GEOM is a valid geometry
%   struct, HMAX and HMAX_E are arrays or scalars indicating the maximum grid size for
%   subdomains and edges, respectively. FILE_NAME is optional and enables output to
%   file (writes a poly file). FID_LOG specifies a message log file handle (negative
%   for gui output or empty for no output).

% Copyright 2013-2024 Precise Simulation, Ltd.
