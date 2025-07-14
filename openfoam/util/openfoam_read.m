%OPENFOAM_READ Read OpenFOAM dictionary.
%
%   [ DICT ] = OPENFOAM_READ( FILE_PATH, CONVERT )
%   Read OpenFOAM dictionary FILE_PATH (in ASCII format).
%
%   The optional CONVERT argument toggles conversion of entries to
%   numeric values (default true).
%
%   Returns a struct DICT with (nested) fields where list entries are
%   prefixed with "L_", and non-valid field names are transformed so
%   that non valid characters are replaced with their hexadecimal
%   value prefixed with "x". For example '#' becomnes "x23". And a
%   variable field name starting with a non-letter char is prefixed
%   with a null characted "x00".
%
%   See also OPENFOAM_WRITE

% Copyright 2013-2025 Precise Simulation, Ltd.
