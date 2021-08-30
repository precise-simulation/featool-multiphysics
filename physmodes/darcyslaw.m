function [ varargout ] = darcyslaw( varargin )
%DARCYSLAW Darcy's law physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help darcyslaw, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'darcyslaw', varargin{:} );
if( ~nargout ), clear varargout; end
