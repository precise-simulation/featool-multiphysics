function [ varargout ] = gobj_block( varargin )
%GOBJ_BLOCK Create block geometry object.
%
%   [ GOBJ ] = GOBJ_BLOCK( XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, TAG, T )
%   Creates a block geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       xmin        scalar  {0}               Minimum x-coordinate
%       xmax        scalar  {1}               Maximum x-coordinate
%       ymin        scalar  {0}               Minimum y-coordinate
%       ymax        scalar  {1}               Maximum y-coordinate
%       zmin        scalar  {0}               Minimum z-coordinate
%       zmax        scalar  {1}               Maximum z-coordinate
%       tag         string  {B1}              Geometry object tag/name
%       t           logical {false}           Triangulate boundary segments

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_block, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_block', varargin{:} );
if( ~nargout ), clear varargout; end
