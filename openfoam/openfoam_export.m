%OPENFOAM_EXPORT Export OpenFOAM grid and problem data.
%
%   [ DICT ] = OPENFOAM_EXPORT( FEA, OPT/VARARGIN ) Exports an
%   OpenFOAM problem defined by the FEA problem struct and options
%   OPT (to the folder specified by opt.meta.casedir, default
%   tempdir/foam_case). Property/value pairs VARARGIN can also be
%   given instead of an OPT struct.
%
%   Optionally returns a struct DICT of OpenFOAM dictionaries.
%
%   See also OPENFOAM, OPENFOAM_DATA, OPENFOAM_WRITE

% Copyright 2013-2025 Precise Simulation, Ltd.
