%ELECTROSTATICS Electrostatics physics mode.
%
%   ELECTROSTATICS Define and generate a electrostatics
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default electrostatics physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @electrostatics );
%
%      2) Add electrostatics physics mode with named dependent variable.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @electrostatics, {'V1'} );
%
%
%   [ PHYS ] = ELECTROSTATICS( SDIM, DVAR, TAG, ISAXI )
%   Optional arguments when used without addphys are a cell array of
%   space dimension variables SDIM (default {'x', 'y'}), dependent
%   variable name DVAR (default {'V'}), mode tag TAG (default
%   'es'), and boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
