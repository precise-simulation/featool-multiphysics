%GRIDREFINE Uniform refinement of a grid.
%
%   [ SOUT ] = GRIDREFINE( SIN, FID ) Creates a uniform refinement of a grid.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and a, s, and b) fields
%       fid         scalar/{1}             File identifier for output ([]=no output)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Output grid or problem struct
%
%   See also GRIDREFINE1, GRIDREFINE2, GRIDREFINE3

% Copyright 2013-2024 Precise Simulation, Ltd.
