function [ varargout ] = planestrain( varargin )
%PLANESTRAIN Plane strain physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help planestrain, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'planestrain', varargin{:} );
if( ~nargout ), clear varargout; end
