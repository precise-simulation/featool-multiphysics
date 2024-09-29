%OPENFOAM_IMPORT Import OpenFOAM solution data.
%
%   [ U, TLIST, VARS/DATA ] = OPENFOAM_IMPORT( FEA / VARARGIN )
%   Imports OpenFOAM solution data from the folder CASEDIR (from
%   numeric time subfolders for dictionaries with internalField).
%
%   If the FEATool problem struct argument FEA is given, the output
%   DATA is converted to compatible [U, TLIST, VARS]. Otherwise the a
%   struct is returned with fields TLIST for the time steps and
%   variable names, DATA with the internalField data (size n_zones x
%   n_timesteps), and ZONES for subdomain names (empty for main domain).
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       applybc     logical {true}               Re-apply (enforce) Dirichlet boundary conditions
%       casedir     string {pwd}                 OpenFOAM case directory
%       interp      scalar {1}                   Interpolation method (cell to vertex)
%                                                  0: no weighting 1: cell volume weighting
%       iswarn      logical {true}               Show warnings
%       tskip       array {[]}                   List of times/directories to skip
%       varmap      n x 2 cell {default}         OpenFOAM to FEATool variable mapping
%       zonemap     n x 2 cell {default}         Subdomain to zone name mapping (zone1, ... zonen)
%
%   See also OPENFOAM_IMPORT, OPENFOAM_READ

% Copyright 2013-2024 Precise Simulation, Ltd.
