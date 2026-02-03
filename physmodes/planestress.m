%PLANESTRESS Plane stress physics mode.
%
%   PLANESTRESS Define and generate a plane stress
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default plane stress physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @planestress );
%
%      2) Add plane stress physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_rectangle};
%      fea = addphys( fea, @planestress, {'u1', 'u2'} );
%
%
%   [ PHYS ] = PLANESTRESS( SDIM, DVAR, TAG ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x', 'y'}), dependent variable names
%   DVAR (default {'u', 'v'}), mode tag TAG (default 'pss').
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
