function [ varargout ] = geom_add_gobj( varargin )
%GEOM_ADD_GOBJ Add geometry object(s) to geom or fea struct.
%
%   [ FEA ] = GEOM_ADD_GOBJ( FEA, GOBJ ) Adds geometry object(s) GOBJ
%   to FEA struct ensuring unique name tags.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_add_gobj, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_add_gobj', varargin{:} );
if( ~nargout ), clear varargout; end
