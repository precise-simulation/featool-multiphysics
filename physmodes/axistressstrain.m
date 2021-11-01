function [ varargout ] = axistressstrain( varargin )
%AXISTRESSSTRAIN Axisymmetric stress-strain physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help axistressstrain, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'axistressstrain', varargin{:} );
if( ~nargout ), clear varargout; end
