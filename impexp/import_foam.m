%IMPORT_FOAM Import OpenFOAM field from file.
%
%   [ DATA ] = IMPORT_FOAM( FILE, FIELD, FID_LOG )
%
%   Imports a numeric OpenFOAM FIELD (default internalField) from FILE
%   (ASCII format). Returns a row array DATA of size n (cells) x m
%   (vector).
%
%   FID_LOG is an optional log file handle for message output
%   (negative for GUI output or empty for no output).
%
%   See also OPENFOAM_IMPORT

% Copyright 2013-2026 Precise Simulation, Ltd.
