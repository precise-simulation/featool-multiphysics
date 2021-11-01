function [ varargout ] = magnetostatics( varargin )
%MAGNETOSTATICS Magnetostatics physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help magnetostatics, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'magnetostatics', varargin{:} );
if( ~nargout ), clear varargout; end
