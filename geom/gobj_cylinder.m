%GOBJ_CYLINDER Create cylinder geometry object.
%
%   [ GOBJ ] = GOBJ_CYLINDER( P, R, L, AX, TAG ) Creates a cylinder
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0,0,0]}         Coordinates of base center point
%       r           scalar  {1}               Cylinder radius
%       l           scalar  {1}               Cylinder length
%       ax          scalar/array {[1,0,0]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       tag         string  {C1}              Geometry object tag/name

% Copyright 2013-2022 Precise Simulation, Ltd.
