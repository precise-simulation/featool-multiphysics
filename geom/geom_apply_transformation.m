function [ varargout ] = geom_apply_transformation( varargin )
%GEOM_APPLY_TRANSFORMATION Apply transformation to geometry objects.
%
%   [ SOUT, NEW_TAGS ] = GEOM_APPLY_TRANSFORMATION( SIN, TAGS, DP, S, TH, AX, P0, DEG )
%   Applies transformations to the geometry objects in finite element
%   or geometry struct SIN with the specified TAGS.
%
%   DP is a (1 x n_sdim) vector specifying linear translation in
%   the coordinate directions (default [0,0,0], no translation).
%
%   S applies a scaling (1 x n_sdim) in the coordinate directions.
%   (default [1,1,1], no scaling)
%
%   TH specifies a rotation in radians around the point P0 (default
%   TH = 0, and P0 = [0,0(,0)]). For 3D geometry objects AX is a vector
%   specifying the axis of rotation (default x-axis [1,0,0]). DEG
%   optionally indicates that TH is given in degrees.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_apply_transformation, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_apply_transformation', varargin{:} );
if( ~nargout ), clear varargout; end
