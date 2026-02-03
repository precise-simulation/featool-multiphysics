%GEOM_ANALYZE Analyze geometry objects.
%
%   [ SOUT ] = GEOM_ANALYZE( SIN, TOL, IS_WARN, IS_SPLIT ) Analyzes
%   geometry (objects) SIN and returns a decomposed geometry with
%   non-overlapping geometry objects. TOL (default 1e-5) specifies the
%   tolerance for CSG boundary operations. The IS_WARN flag (default
%   false) enables throwing warnings if the geometry objects cannot be
%   combined. IS_SPLIT (default true) returns separate objects, or a
%   single combined one (false, only with geomtool).

% Copyright 2013-2026 Precise Simulation, Ltd.
