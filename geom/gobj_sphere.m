function [ varargout ] = gobj_sphere( varargin )
%GOBJ_SPHERE Create sphere geometry object.
%
%   [ GOBJ ] = GOBJ_SPHERE( P, R, AX, TAG, N_S, T ) Creates a sphere
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       r           scalar  {1}               Sphere radius
%       ax          scalar/array {1}          Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       tag         string  {S1}              Geometry object tag/name
%       t           logical {false}           Triangulate boundary segments

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_sphere, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_sphere', varargin{:} );
if( ~nargout ), clear varargout; end
