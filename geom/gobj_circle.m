function [ varargout ] = gobj_circle( varargin )
%GOBJ_CIRCLE Create circle geometry object.
%
%   [ GOBJ ] = GOBJ_CIRCLE( P, R, TAG ) Creates a circle
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0]}            Coordinates of center point
%       r           scalar {1}                Circle radius
%       tag         string {C1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_circle, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_circle', varargin{:} );
if( ~nargout ), clear varargout; end
