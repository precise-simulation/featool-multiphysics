%PARSEPHYS Check and parse physics modes.
%
%   [ SOUT ] = PARSEPHYS( SIN, F_ADD ) Checks and parses the physics modes
%   in the given struct SIN.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Problem definition or physics mode struct
%       f_add       logical                Flag to expand and add physics mode
%       f_constr    logical                Flag to add (ns/br) pressure constraint
%                                          data to problem struct fields (default true)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Parsed problem or physics mode struct
%
%   See also PARSEPROB, PARSEPROB0

% Copyright 2013-2024 Precise Simulation, Ltd.
