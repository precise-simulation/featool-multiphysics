%FEACACHE FEATool data cache.
%
%   [ OUT ] = FEACACHE( VAR, SLOT, VAL ) Cache and retrieve data. VAR
%   designates a cache variable name, an integer SLOT retrieves, or
%   puts the argument VAL to that slot (cache{VAR,SLOT} = VAL).
%
%   Output OUT is either a requested VAR/SLOT pair, or with not input
%   arguments a struct ST_CACHE of cached expressions used by the
%   feacache, and evalexpr0, parseexpr, subcoef, and symvarbe.
%
%   FEACACHE(0) Clears the cache.

% Copyright 2013-2022 Precise Simulation, Ltd.
