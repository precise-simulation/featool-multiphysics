%GOBJ_CONE Create cone geometry object.
%
%   [ GOBJ ] = GOBJ_CONE( P, R1, R2, L, AX, TAG ) Creates a cone
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0,0,0]}         Coordinates of base center point
%       r1          scalar  {1}               Base cone radius
%       r2          scalar  {0.5}             End cone radius
%       l           scalar  {1}               Cylinder length
%       ax          scalar/array {[1,0,0]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       tag         string  {C1}              Geometry object tag/name

% Copyright 2013-2024 Precise Simulation, Ltd.
