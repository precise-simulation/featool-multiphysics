%LINEARELASTICITY Linear elasticity physics mode.
%
%   LINEARELASTICITY Define and generate a solid stress-strain
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default solid stress-strain physics mode.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_cylinder};
%      fea = addphys( fea, @linearelasticity );
%
%      2) Add solid stress-strain physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @linearelasticity, {'u1', 'u2', 'u3'} );
%
%
%   [ PHYS ] = LINEARELASTICITY( SDIM, DVAR, TAG ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x', 'y', 'z'}), dependent variable names
%   DVAR (default {'u', 'v', 'w'}), mode tag TAG (default 'el').
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
