%MAGNETOSTATICS Magnetostatics physics mode.
%
%   MAGNETOSTATICS Define and generate a magnetostatics
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default magnetostatics physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @magnetostatics );
%
%      2) Add magnetostatics physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @magnetostatics, {'A1', 'A2', 'A3'} );
%
%
%   [ PHYS ] = MAGNETOSTATICS( SDIM, DVAR, TAG, ISAXI ) Optional
%   arguments when used without addphys are a cell array of space
%   dimension variables C_SDIM (default {'x', 'y'}), dependent
%   variable name DVAR (default {'Vm'} in 2D, {'Ax', 'Ay', 'Az'} in
%   3D), mode tag TAG (default 'ms'), and boolean ISAXI to specify 2D
%   axisymmetric formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
