%GRIDROTATE Rotate grid.
%
%   [ GRID ] = GRIDROTATE( GRID, TH, AX ) Applies rotation angle TH
%   (radians) to grid points around the origin (in 3D around axis AX:
%   1 = x-axis (default), 2 = y-axis, 3 = z-axis). The input argument
%   GRID can either be a grid struct or simply an array of
%   coordinates. Normals will also be recalculated in the case of a
%   grid struct with boundary information.
%
%   See also GRIDEXTRUDE, GRIDMERGE, GRIDREVOLVE, GRIDSCALE

% Copyright 2013-2022 Precise Simulation, Ltd.
