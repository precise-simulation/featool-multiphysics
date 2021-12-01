function [ varargout ] = poisson( varargin )
%POISSON Poisson equation physics mode definition.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help poisson, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'poisson', varargin{:} );
if( ~nargout ), clear varargout; end
