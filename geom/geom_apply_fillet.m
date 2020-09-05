function [ varargout ] = geom_apply_fillet( varargin )
%GEOM_APPLY_FILLET Apply fillet to geometry objects.
%
%   [ SOUT, NEW_TAGS ] = GEOM_APPLY_FILLET( SIN, TAGS, R ) Applies
%   fillets to the geometry objects in finite element or geometry
%   struct SIN with the specified TAGS. R specifies the fillet radius
%   to apply.

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_apply_fillet, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_apply_fillet', varargin{:} );
if( ~nargout ), clear varargout; end
