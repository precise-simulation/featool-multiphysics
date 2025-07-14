%DARCYSLAW Darcy's law physics mode.
%
%   DARCYSLAW Define and generate a Darcy's law physics mode struct.
%   Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default Darcy's law physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @darcyslaw );
%
%      2) Add Darcy's law physics mode with named dependent variable.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @darcyslaw, {'p0'} );
%
%
%   [ PHYS ] = DARCYSLAW( SDIM, DVAR, TAG, ISAXI ) Optional
%   arguments when used without addphys are a cell array of space
%   dimension variables SDIM (default {'x', 'y'}), dependent
%   variable name DVAR (default {'p'}), mode tag TAG (default
%   'dl'), and boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2025 Precise Simulation, Ltd.
