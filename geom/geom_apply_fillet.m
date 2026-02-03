%GEOM_APPLY_FILLET Apply fillet to geometry objects.
%
%   [ SOUT, NEW_TAGS ] = GEOM_APPLY_FILLET( SIN, TAGS, R, IND_F )
%   Applies fillets to the geometry objects in finite element or
%   geometry struct SIN with the specified TAGS. R specifies the
%   fillet radius to apply, and IND_F optionally selects faces to
%   apply fillets to (default all). Negative IND_F values indicates
%   vertices/edges in 2D/3D respectively.

% Copyright 2013-2026 Precise Simulation, Ltd.
