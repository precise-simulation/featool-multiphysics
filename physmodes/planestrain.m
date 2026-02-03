%PLANESTRAIN Plane strain physics mode.
%
%   PLANESTRAIN Define and generate a plane strain
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default plane strain physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @planestrain );
%
%      2) Add plane strain physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_rectangle};
%      fea = addphys( fea, @planestrain, {'u1', 'u2'} );
%
%
%   [ PHYS ] = PLANESTRAIN( SDIM, DVAR, TAG ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x', 'y'}), dependent variable names
%   DVAR (default {'u', 'v'}), mode tag TAG (default 'psn').
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
