%GOBJ_TORUS Create torus geometry object.
%
%   [ GOBJ ] = GOBJ_TORUS( P, RC, RT, AX, TAG ) Creates an torus
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       rc          scalar  {1}               Torus radius (tube center to torus center)
%       rt          scalar  {0.25}            Torus tube radius
%       ax          scalar/array {[1,0,0]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,1])
%       tag         string  {T1}              Geometry object tag/name

% Copyright 2013-2024 Precise Simulation, Ltd.
