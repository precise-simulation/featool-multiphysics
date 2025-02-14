%GEOM_APPLY_FORMULA Apply formula to a geometry struct.
%
%   [ SOUT, STAT ] = GEOM_APPLY_FORMULA( SIN, FORMULA, IS_WARN )
%   Applies a FORMULA to geometry objects in the fea or geometry
%   struct SIN. The FORMULA is a string containing combinations of
%   geometry object tags and operations ('+' union/join, '-' subtract,
%   '&' intersect). The status return code STAT is zero if the
%   operation has been applied sucessfully, or a positive integer
%   corresponding to the number of failed geometry object
%   combinations. The IS_WARN flag (default false) enables throwing
%   warnings if the geometry objects cannot be combined.

% Copyright 2013-2025 Precise Simulation, Ltd.
