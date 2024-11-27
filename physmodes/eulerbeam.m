%EULERBEAM 1D Euler-Bernoulli beam physics mode.
%
%   EULERBEAM Define and generate a Euler beam
%   physics mode struct. Used with the ADDPHYS function.
%
%   Example:
%
%      1) Add default eulerbeam physics mode.
%
%      fea.sdim = {'x'};
%      fea.geom.objects = {gobj_line};
%      fea = addphys( fea, @eulerbeam );
%
%
%   [ PHYS ] = EULERBEAM( SDIM, DVAR, TAG ) Optional arguments
%   when used without addphys are a cell array of space dimension
%   variables SDIM (default {'x'}), dependent variable name DVAR
%   (default {'v'}), mode tag TAG (default 'eb').
%
%   See also ADDPHYS

% Copyright 2013-2024 Precise Simulation, Ltd.
