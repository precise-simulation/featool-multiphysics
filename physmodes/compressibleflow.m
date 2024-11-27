%COMPRESSIBLEFLOW Compressible viscous flow physics mode.
%
%   COMPRESSIBLEFLOW Define and generate a compressible flow
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default compressible flow physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @compressibleflow );
%
%      2) Add compressible flow physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @compressibleflow, {'T0', 'u2', 'u3', 'p'} );
%
%
%   [ PHYS ] = COMPRESSIBLEFLOW( SDIM, DVAR, TAG, ISAXI )
%   Optional arguments when used without addphys are a cell array of
%   space dimension variables SDIM (default {'x', 'y'}), dependent
%   variable names DVAR (default {'T', 'u', 'v', 'p'}), mode tag
%   TAG (default 'cf'), and boolean ISAXI to specify 2D axisymmetric
%   formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
