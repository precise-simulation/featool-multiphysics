%DELCELLS Delete cells from grid.
%
%   [ SOUT ] = DELCELLS( SIN, SELECTION ) Deletes cells from a grid.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and s) fields
%       selection   array                  Array with indicies or string expression
%                                          indicating cells to delete
%       bdr         logical/{true}         Reconstruct boundary field
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Output grid or problem struct

% Copyright 2013-2025 Precise Simulation, Ltd.
