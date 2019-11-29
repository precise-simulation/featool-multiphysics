function [ varargout ] = openfoam( varargin )
%OPENFOAM Exports, solves, and imports an OpenFOAM CFD problem.
%
%   [ U, TLIST, VARS ] = OPENFOAM( PROB, VARARGIN ) Export, solves, or
%   imports the solved problem described in the PROB finite element
%   struct using the OpenFOAM CFD solver. Accepts the following
%   property/value pairs, moreover also accepts the openfoam_data
%   property/value pairs to set the OpenFOAM dicts during export.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       mode        check, export, solve, import Command mode(s) to call (default all)
%       data        default                      Default OpenFOAM data and parameter dict
%       init        scalar {[]}                  Initial values   []: init expressions
%                                                  i: use solution -1: potential flow
%       turb        struct                       Turbulence data struct fields
%                                                  model, inlet [k,e/o], wallfcn [1/0]
%       interp      scalar {2}                   Interpolate solution to grid points
%                                                  1: no weighting 2: cell volume weighting
%       control     logical {0}                  Show solver control panel.
%       casedir     default                      OpenFOAM case directory
%       foamdir     default                      OpenFOAM installation directory
%       logfname    default                      OpenFOAM log/output filename
%       fid/logfid  scalar {1}                   Log file/message output file handle
%       hax         handle                       Axis handle to plot convergence
%       naxts       scalar {250}                 Maximum number of time steps to plot
%
%   See also OPENFOAM_DATA

% Copyright 2013-2019 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help openfoam, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'openfoam', varargin{:} );
if( ~nargout ), clear varargout; end
