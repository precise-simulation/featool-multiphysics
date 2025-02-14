%HEATTRANSFER Heat transfer physics mode definition.
%
%   HEATTRANSFER Define and generate a heat transfer mode struct.
%   Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default heat transfer physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @heattransfer );
%
%      2) Add heat transfer physics mode with named dependent variable.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @heattransfer, {'T2'} );
%
%
%   [ PHYS ] = HEATTRANSFER( SDIM, DVAR, TAG, ISAXI ) Optional
%   arguments when used without addphys are a cell array of space
%   dimension variables SDIM (default {'x', 'y'}), dependent variable
%   name DVAR (default {'T'}), mode tag TAG (default 'ht'), and
%   boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2025 Precise Simulation, Ltd.
