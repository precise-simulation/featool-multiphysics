function [ varargout ] = gobj_cylinder( varargin )
%GOBJ_CYLINDER Create cylinder geometry object.
%
%   [ GOBJ ] = GOBJ_CYLINDER( P, R, L, AX, TAG, N_S, T ) Creates a cylinder
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0,0,0]}         Coordinates of base center point
%       r           scalar  {1}               Cylinder radius
%       l           scalar  {1}               Cylinder length
%       ax          scalar/array {[1,0,0]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       tag         string  {C1}              Geometry object tag/name
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       t           logical {false}           Triangulate boundary segments

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_cylinder, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_cylinder', varargin{:} );
if( ~nargout ), clear varargout; end
