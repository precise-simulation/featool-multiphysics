function [ varargout ] = geom_add_gobj( varargin )
%GEOM_ADD_GOBJ Add geometry object(s) to geom or fea struct.

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_add_gobj, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_add_gobj', varargin{:} );
if( ~nargout ), clear varargout; end
