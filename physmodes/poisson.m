%POISSON Poisson equation physics mode definition.
%
%   POISSON Define and generate a Poisson physics mode struct.
%   Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default Poisson physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @poisson );
%
%      2) Add Poisson physics mode with named dependent variable.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @poisson, {'U'} );
%
%
%   [ PHYS ] = POISSON( SDIM, DVAR, TAG, ISAXI ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x', 'y'}), dependent variable name DVAR
%   (default {'u'}), mode tag TAG (default 'poi'), and boolean ISAXI
%   to specify 2D axisymmetric formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
