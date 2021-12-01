function [ varargout ] = gobj_ellipsoid( varargin )
%GOBJ_ELLIPSOID Create ellipsoid geometry object.
%
%   [ GOBJ ] = GOBJ_ELLIPSOID( P, RX, RY, RZ, AX, TAG, N_S, T )
%   Creates an ellipsoid geometry object. Accepts the following input
%   parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       rx          scalar  {1}               Ellipsoid x-radius
%       ry          scalar  {2}               Ellipsoid y-radius
%       rz          scalar  {3}               Ellipsoid z-radius
%       ax          scalar/array {[0,0,1]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis rotation vector (ex. [1,1,1])
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       tag         string  {E1}              Geometry object tag/name
%       t           logical {false}           Triangulate boundary segments

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_ellipsoid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_ellipsoid', varargin{:} );
if( ~nargout ), clear varargout; end
