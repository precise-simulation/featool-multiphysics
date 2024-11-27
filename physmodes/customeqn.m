%CUSTOMEQN Custom equation physics mode definition.
%
%   CUSTOMEQN Define and generate a custom physics mode struct.
%   Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default custom equation physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @customeqn );
%
%      2) Add custom equation physics mode with 3 named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @customeqn, {'u1', 'u2', 'u3'} );
%
%
%   [ PHYS ] = CUSTOMEQN( SDIM, DVAR, TAG, ISAXI ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x', 'y'}), dependent variable names DVAR
%   (default {'u'}), mode tag TAG (default 'ce'), and boolean ISAXI
%   to specify 2D axisymmetric formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
