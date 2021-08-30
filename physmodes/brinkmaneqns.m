function [ varargout ] = brinkmaneqns( varargin )
%BRINKMANEQNS Brinkman equations physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help brinkmaneqns, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'brinkmaneqns', varargin{:} );
if( ~nargout ), clear varargout; end
