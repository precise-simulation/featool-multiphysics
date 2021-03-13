function [ varargout ] = compressibleeuler( varargin )
%COMPRESSIBLEEULER Compressible Euler equations physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help compressibleeuler, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'compressibleeuler', varargin{:} );
if( ~nargout ), clear varargout; end
