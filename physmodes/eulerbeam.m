function [ varargout ] = eulerbeam( varargin )
%EULERBEAM 1D Euler-Bernoulli beam physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help eulerbeam, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'eulerbeam', varargin{:} );
if( ~nargout ), clear varargout; end
