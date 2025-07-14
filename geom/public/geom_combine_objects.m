%GEOM_COMBINE_OBJECTS Combine geometry objects.
%
%   [ GEOM, NEW_TAG, STAT ] = GEOM_COMBINE_OBJECTS( GEOM, TAG1, TAG2, OP, IS_WARN )
%   Combines geometry objects with TAG1 and TAG2 through application
%   of OP, where OP is a string character indicating operation to
%   perform ('+' union/join, '-' subtract, '&' intersect). The status
%   return code STAT is zero if the operation has been applied
%   sucessfully, and a positive integer to indicate error. The IS_WARN
%   flag (default false) enables throwing warnings if the geometry
%   objects cannot be combined.

% Copyright 2013-2025 Precise Simulation, Ltd.
