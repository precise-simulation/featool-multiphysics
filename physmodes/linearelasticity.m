function [ varargout ] = linearelasticity( varargin )
%LINEARELASTICITY Linear elasticity physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help linearelasticity, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'linearelasticity', varargin{:} );
if( ~nargout ), clear varargout; end
