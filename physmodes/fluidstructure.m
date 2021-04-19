function [ varargout ] = fluidstructure( varargin )
%FLUIDSTRUCTURE Fluid-structure interaction physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fluidstructure, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fluidstructure', varargin{:} );
if( ~nargout ), clear varargout; end
