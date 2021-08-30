function [ varargout ] = geom_apply_chamfer( varargin )
%GEOM_APPLY_FILLET Apply chamfer to geometry objects.
%
%   [ SOUT, NEW_TAGS ] = GEOM_APPLY_CHAMFER( SIN, TAGS, D, IND_F )
%   Applies chamfers to the geometry objects in finite element or
%   geometry struct SIN with the specified TAGS. D specifies the
%   chamfer distance to apply, and IND_F optionally selects faces to
%   apply chamfers to (default all).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_apply_chamfer, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_apply_chamfer', varargin{:} );
if( ~nargout ), clear varargout; end
