%FLUIDSTRUCTURE Fluid-structure interaction physics mode.
%
%   FLUIDSTRUCTURE Define and generate a fluid-structure interaction
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default fluid-structure interaction physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_rectangle};
%      fea = addphys( fea, @fluidstructure );
%
%
%   [ PHYS ] = FLUIDSTRUCTURE( SDIM, DVAR, TAG, ISAXI )
%   Optional arguments when used without addphys are a cell array of
%   space dimension variables SDIM (default {'x', 'y'}), dependent
%   variable names DVAR (default {'u', 'v', 'dx', 'dy'}), mode tag
%   TAG (default 'fsi'), and boolean ISAXI to specify 2D axisymmetric
%   formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
