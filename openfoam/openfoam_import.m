%OPENFOAM_IMPORT Import OpenFOAM solution data.
%
%   [ U, TLIST, VARS/DATA ] = OPENFOAM_IMPORT( FEA / VARARGIN )
%   Imports OpenFOAM solution data from the folder CASEDIR
%   (from numeric time subfolders for dictionaries with internalField).
%
%   If the FEATool problem struct argument FEA is given, the output
%   DATA is converted to compatible [U, TLIST, VARS]. Otherwise a
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
%       varmap      N x 2 cell {default}         OpenFOAM to FEATool variable mapping
%       zonemap     N x 2 cell {default}         Subdomain to zone name mapping (zone1, ... zoneN)
%
%   Examples:
%
%      1) Laminar steady Hagen-Poiseuille flow in a channel.
%
%      fea = ex_navierstokes1( 'iplot', false );
%      casedir = 'C:\temp\oftest';
%      openfoam_export( fea, 'casedir', casedir )
%      ls(casedir)
%      data = openfoam_import( fea, 'casedir', casedir )
%
%   See also OPENFOAM_EXPORT, OPENFOAM_READ, EX_NAVIERSTOKES1

% Copyright 2013-2025 Precise Simulation, Ltd.
