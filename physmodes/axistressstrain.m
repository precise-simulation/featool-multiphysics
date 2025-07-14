%AXISTRESSSTRAIN Axisymmetric stress-strain physics mode.
%
%   AXISTRESSSTRAIN Define and generate a axisymmetric stress-strain
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default axisymmetric stress-strain physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @axistressstrain );
%
%      2) Add axisymmetric stress-strain physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_rectangle};
%      fea = addphys( fea, @axistressstrain, {'ur', 'uz'} );
%
%
%   [ PHYS ] = AXISTRESSSTRAIN( SDIM, DVAR, TAG ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'r', 'z'}), dependent variable name
%   DVAR (default {'u', 'w'}), mode tag TAG (default 'css').
%
%   See also ADDPHYS

% Copyright 2013-2025 Precise Simulation, Ltd.
