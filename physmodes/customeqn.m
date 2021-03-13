function [ varargout ] = customeqn( varargin )
%CUSTOMEQN Custom equation physics mode definition.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help customeqn, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'customeqn', varargin{:} );
if( ~nargout ), clear varargout; end
