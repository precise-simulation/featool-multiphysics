function [ varargout ] = impexp_foam( varargin )
%IMPEXP_FOAM Import/export OpenFOAM format.
%
%   [ U, VARS ] = IMPEXP_FOAM( CASEDIR, MODE, PROB, DICT, FID_LOG )
%   Import and export of OpenFOAM ASCII mesh and problem data.
%   Returns the solution vector U (after import) and optional
%   turbulence quantitiy variables in the VARS struct.
%
%   CASEDIR specifies the location of OpenFOAM case files.
%
%   MODE can either be a string indicating IMPORT (default) or EXPORT.
%
%   PROB is a valid FEATool FEA problem struct.
%
%   DICT are the OpenFOAM dict structs to write.
%
%   FID_LOG is an optional log file handle for message output
%   (negative for GUI output or empty for no output).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_foam, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_foam', varargin{:} );
if( ~nargout ), clear varargout; end
