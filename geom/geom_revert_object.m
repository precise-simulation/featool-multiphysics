function [ varargout ] = geom_revert_object( varargin )
%GEOM_REVERT_OBJECT Revert/undo geometry object operations.
%
%   [ GEOM ] = GEOM_REVERT_OBJECT( GEOM, TAG ) Reverts/undos geometry
%   object operations by replacing the geometry object TAG with the
%   objects in its children field (with uniqified tags).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_revert_object, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_revert_object', varargin{:} );
if( ~nargout ), clear varargout; end
