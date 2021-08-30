function [ varargout ] = gobj_ellipse( varargin )
%GOBJ_ELLIPSE Create ellipse geometry object.
%
%   [ GOBJ ] = GOBJ_ELLIPSE( P, RX, RY, TAG ) Creates an ellipse
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0]}            Coordinates of center point
%       rx          scalar {1}                Radius along x-axis
%       ry          scalar {0.5}              Radius along y-axis
%       tag         string {E1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_ellipse, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_ellipse', varargin{:} );
if( ~nargout ), clear varargout; end
