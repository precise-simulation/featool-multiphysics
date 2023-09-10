%RMPHYS Remove physics modes from problem struct.
%
%   [ PROB ] = RMPHYS( PROB, SEL ) Removes physics modes from a
%   problem struct PROB. SEL is either an array of integers or cell
%   array of physics mode tags, to indicate modes to remove.
%   Calls parsephys and parseprob to reset dvar, sfun, eqn, and bdr
%   fields.
%
%   See also ADDPHYS

% Copyright 2013-2023 Precise Simulation, Ltd.
