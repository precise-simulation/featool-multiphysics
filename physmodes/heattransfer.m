function [ varargout ] = heattransfer( varargin )
%HEATTRANSFER Heat transfer physics mode definition.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help heattransfer, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'heattransfer', varargin{:} );
if( ~nargout ), clear varargout; end
