function [ varargout ] = conductivemediadc( varargin )
%CONDUCTIVEMEDIADC Conductive media DC physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help conductivemediadc, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'conductivemediadc', varargin{:} );
if( ~nargout ), clear varargout; end
