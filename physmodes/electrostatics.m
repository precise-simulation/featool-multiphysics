function [ varargout ] = electrostatics( varargin )
%ELECTROSTATICS Electrostatics physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help electrostatics, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'electrostatics', varargin{:} );
if( ~nargout ), clear varargout; end
