function [ varargout ] = convectiondiffusion( varargin )
%CONVECTIONDIFFUSION Convection and diffusion physics mode definition.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help convectiondiffusion, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'convectiondiffusion', varargin{:} );
if( ~nargout ), clear varargout; end
