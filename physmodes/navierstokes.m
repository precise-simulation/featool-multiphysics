function [ varargout ] = navierstokes( varargin )
%NAVIERSTOKES Navier-Stokes equations physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help navierstokes, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'navierstokes', varargin{:} );
if( ~nargout ), clear varargout; end
