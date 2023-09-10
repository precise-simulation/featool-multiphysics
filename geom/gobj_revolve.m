%GOBJ_REVOLVE Revolve face/geometry object.
%
%   [ GOBJ, STAT ] = GOBJ_REVOLVE( GOBJ, TH, P, V, PIT, FACE_WPL )
%   Revolve geometry object or face from GOBJ a an angle TH (degrees)
%   around axis defined by point P and direction vector V. PIT is
%   optionally the pitch (distance per revolution) to offset the
%   revolution in the axial direction (to form a coil/helix).
%   FACE_WPL can in 3D indicate a face to extrude or in 2D is a
%   workplane struct defined by (p, n, t).

% Copyright 2013-2023 Precise Simulation, Ltd.
