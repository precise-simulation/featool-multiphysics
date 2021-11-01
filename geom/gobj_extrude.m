function [ varargout ] = gobj_extrude( varargin )
%GOBJ_EXTRUDE Extrude face/geometry object.
%
%   [ GOBJ, STAT ] = GOBJ_EXTRUDE( GOBJ, D, V, FACE_WPL )
%   Extrude geometry object or face from GOBJ a distance D in the
%   direction of the vector V. FACE_WPL can in 3D indicate a face to
%   extrude or in 2D is a workplane struct defined by (p, n, t).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_extrude, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_extrude', varargin{:} );
if( ~nargout ), clear varargout; end
