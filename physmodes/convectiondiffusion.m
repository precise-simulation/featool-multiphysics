%CONVECTIONDIFFUSION Convection and diffusion physics mode definition.
%
%   CONVECTIONDIFFUSION Define and generate a convection and diffusion
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default convection and diffusion physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @convectiondiffusion );
%
%      2) Add convection and diffusion physics mode with 3 named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @convectiondiffusion, {'c1', 'c2', 'c3'} );
%
%
%   [ PHYS ] = CONVECTIONDIFFUSION( SDIM, DVAR, TAG, ISAXI ) Optional
%   arguments when used without addphys are a cell array of space
%   dimension variables SDIM (default {'x', 'y'}), dependent variable
%   names DVAR (default {'c'}), mode tag TAG (default 'cd'), and
%   boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
