%CONDUCTIVEMEDIADC Conductive media DC physics mode.
%
%   CONDUCTIVEMEDIADC Define and generate a conductive media DC
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default conductive media DC physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @conductivemediadc );
%
%      2) Add conductive media DC physics mode with named dependent variable.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @conductivemediadc, {'V1'} );
%
%
%   [ PHYS ] = CONDUCTIVEMEDIADC( SDIM, DVAR, TAG, ISAXI )
%   Optional arguments when used without addphys are a cell array of
%   space dimension variables SDIM (default {'x', 'y'}), dependent
%   variable name DVAR (default {'V'}), mode tag TAG (default
%   'dc'), and boolean ISAXI to specify 2D axisymmetric formulation (default
%   false).
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
