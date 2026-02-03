%GOBJ_ELLIPSOID Create ellipsoid geometry object.
%
%   [ GOBJ ] = GOBJ_ELLIPSOID( P, RX, RY, RZ, AX, TAG )
%   Creates an ellipsoid geometry object. Accepts the following
%   input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       rx          scalar  {1}               Ellipsoid x-radius
%       ry          scalar  {2}               Ellipsoid y-radius
%       rz          scalar  {3}               Ellipsoid z-radius
%       ax          scalar/array {[0,0,1]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis rotation vector (ex. [1,1,1])
%       tag         string  {E1}              Geometry object tag/name

% Copyright 2013-2026 Precise Simulation, Ltd.
