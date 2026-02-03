%NONNEWTONIAN Non-Newtonian flow physics mode.
%
%   NONNEWTONIAN Define and generate a non-Newtonian flow
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default non-Newtonian flow physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @nonnewtonian );
%
%      2) Add non-Newtonian flow physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @nonnewtonian, {'u1', 'u2', 'u3', 'p'} );
%
%
%   [ PHYS ] = NONNEWTONIAN( SDIM, DVAR, TAG, ISAXI ) Optional
%   arguments when used without addphys are a cell array of space
%   dimension variables SDIM (default {'x', 'y'}), dependent variable
%   names DVAR (default {'u', 'v', 'p'}), mode tag TAG (default 'nn'),
%   and boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
