%COMPRESSIBLEEULER Compressible Euler equations physics mode.
%
%   COMPRESSIBLEEULER Define and generate a compressible Euler equations
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default compressible Euler equations physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @compressibleeuler );
%
%      2) Add compressible Euler equations physics mode with named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @compressibleeuler, {'rho0', 'u2', 'u3', 'p'} );
%
%
%   [ PHYS ] = COMPRESSIBLEEULER( SDIM, DVAR, TAG, ISAXI )
%   Optional arguments when used without addphys are a cell array of
%   space dimension variables SDIM (default {'x', 'y'}), dependent
%   variable names DVAR (default {'rho', 'u', 'v', 'p'}), mode tag
%   TAG (default 'ee'), and boolean ISAXI to specify 2D axisymmetric
%   formulation (default false).
%
%   See also ADDPHYS

% Copyright 2013-2026 Precise Simulation, Ltd.
