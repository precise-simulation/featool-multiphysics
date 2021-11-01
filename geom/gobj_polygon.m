function [ varargout ] = gobj_polygon( varargin )
%GOBJ_POLYGON Create polygon geometry object.
%
%   [ GOBJ ] = GOBJ_POLYGON( P, TAG ) Creates a polygon
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0;1 0;0 1]}    Polygon vertex points
%       tag         string {P1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_polygon, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_polygon', varargin{:} );
if( ~nargout ), clear varargout; end
