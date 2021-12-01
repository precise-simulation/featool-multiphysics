function [ varargout ] = planestress( varargin )
%PLANESTRESS Plane stress physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help planestress, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'planestress', varargin{:} );
if( ~nargout ), clear varargout; end
