function [ varargout ] = geom_apply_transformation( varargin )
%GEOM_APPLY_TRANSFORMATION Apply transformation to geometry objects.
%
%   [ SOUT ] = GEOM_APPLY_TRANSFORMATION( SIN, TAGS, DP, S, TH, AX )
%   Applies transformations to the geometry objects with specified
%   TAGS. DP is a 1 x n_sdim array specifying linear translations in
%   the coordinate directions. S applies a scaling (1 x n_sdim), and
%   TH specifies a rotation in radians around the origin in 2D, and
%   AX-axis (1=x default, 2=y, 3=z) in 3D.

% Copyright 2013-2019 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_apply_transformation, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_apply_transformation', varargin{:} );
if( ~nargout ), clear varargout; end
