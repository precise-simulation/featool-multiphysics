function [ varargout ] = addphys( varargin )
%ADDPHYS Parse and add physics mode to problem struct.
%
%   [ PROB ] = ADDPHYS( PROB, F_PHYS, C_DVAR, F_OUT ) Parses and adds
%   the physics mode with function name or handle F_PHYS to the
%   problem struct PROB. C_DVAR prescribes the given instead of
%   default dependent variables names, and the F_OUT flag adds the
%   mode to the PROB struct (default) or optionally simply returns the
%   computed physics mode.
%
%   See also PARSEPHYS

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help addphys, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'addphys', varargin{:} );
if( ~nargout ), clear varargout; end
