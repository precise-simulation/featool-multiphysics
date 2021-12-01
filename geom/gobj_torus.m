function [ varargout ] = gobj_torus( varargin )
%GOBJ_TORUS Create torus geometry object.
%
%   [ GOBJ ] = GOBJ_TORUS( P, RC, RT, AX, TAG, N_S ) Creates an torus
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       rc          scalar  {1}               Torus radius (tube center to torus center)
%       rt          scalar  {0.25}            Torus tube radius
%       ax          scalar/array {[1,0,0]}    Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,1])
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       tag         string  {T1}              Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_torus, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_torus', varargin{:} );
if( ~nargout ), clear varargout; end
