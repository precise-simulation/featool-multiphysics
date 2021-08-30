function [ varargout ] = gobj_point( varargin )
%GOBJ_POINT Create point geometry object.
%
%   [ GOBJ ] = GOBJ_POINT( P, TAG ) Creates a point
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0]}            Coordinates of point
%       tag         string {P1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_point, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_point', varargin{:} );
if( ~nargout ), clear varargout; end
