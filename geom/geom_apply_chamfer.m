%GEOM_APPLY_CHAMFER Apply chamfer to geometry objects.
%
%   [ SOUT, NEW_TAGS ] = GEOM_APPLY_CHAMFER( SIN, TAGS, D, IND_F )
%   Applies chamfers to the geometry objects in finite element or
%   geometry struct SIN with the specified TAGS. D specifies the
%   chamfer distance to apply, and IND_F optionally selects faces to
%   apply chamfers to (default all). Negative IND_F values indicates
%   vertices/edges in 2D/3D respectively.

% Copyright 2013-2024 Precise Simulation, Ltd.
