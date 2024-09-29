%GOBJ_POLYHEDRON Create polyhedron geometry object.
%
%   [ GOBJ ] = GOBJ_POLYHEDRON( P, TAG ) Creates a polyhedron geometry
%   object by finding the convex hull of the points P. Accepts the
%   following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array [n_p,3]             Polyhedron points (default unit tetrahedron)
%       tag         string  {P1}              Geometry object tag/name
%
%   See also CONVHULL

% Copyright 2013-2024 Precise Simulation, Ltd.
