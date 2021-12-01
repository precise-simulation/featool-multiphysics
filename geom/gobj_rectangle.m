function [ varargout ] = gobj_rectangle( varargin )
%GOBJ_RECTANGLE Create rectangle geometry object.
%
%   [ GOBJ ] = GOBJ_RECTANGLE( XMIN, XMAX, YMIN, YMAX, TAG ) Creates a rectangle
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       xmin        scalar {0}                Minimum x-coordinate
%       xmax        scalar {1}                Maximum x-coordinate
%       ymin        scalar {0}                Minimum y-coordinate
%       ymax        scalar {1}                Maximum y-coordinate
%       tag         string {R1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_rectangle, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_rectangle', varargin{:} );
if( ~nargout ), clear varargout; end
