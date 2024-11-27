%OPENFOAM_SOLVE Call external OpenFOAM solver.
%
%   OPENFOAM_SOLVE( VARARGIN ) Call external OpenFOAM solver
%   from the given case directory and monitor solution.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       control     logical {false}              Show solver control panel.
%       casedir     string {current directory}   OpenFOAM case directory
%       nsdim       integer {3}                  Number of space dimensions (2, 2.5 (axi), or 3)
%       foamdir     default                      OpenFOAM installation directory
%       logfname    default                      OpenFOAM log/output filename
%       fid/logfid  scalar {1}                   Log file/message output file handle
%       hax         handle                       Axis handle to plot convergence
%       pmaxts      integer {500}                Maximum number of time steps to print/plot
%       init        scalar {[]}                  Initial values:   []: init expressions
%                                                  i: use solution -1: potential flow
%       tolres      scalar {1e-4/U, 1e-2/p}/vector  Stopping criteria for residuals
%
%   See also OPENFOAM

% Copyright 2013-2024 Precise Simulation, Ltd.
