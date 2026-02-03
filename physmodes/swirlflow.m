%SWIRLFLOW Axisymmetric swirl flow physics mode.
%
%   SWIRLFLOW Define and generate a swirl flow
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default swirl flow physics mode.
%
%      fea.sdim = {'r', 'z'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @swirlflow );
%
%      2) Add swirl flow physics mode with named dependent variables.
%
%      fea.sdim = {'r', 'z'};
%      fea.geom.objects = {gobj_rectangle};
%      fea = addphys( fea, @swirlflow, {'u1', 'u2', 'u3', 'p'} );
%
%
%   [ PHYS ] = SWIRL( SDIM, DVAR, TAG ) Optional arguments when used
%   without addphys are a cell array of space dimension variables SDIM
%   (default {'r', 'z'}), dependent variable names DVAR (default {'u',
%   'v', 'w', 'p'}), mode tag TAG (default 'sw').
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
