function [ varargout ] = geom_apply_formula( varargin )
%GEOM_APPLY_FORMULA Apply formula to a geometry struct.
%
%   [ SOUT, STAT ] = GEOM_APPLY_FORMULA( SIN, FORMULA, TOL, IS_WARN )
%   Applies a FORMULA to geometry objects in the fea or geometry
%   struct SIN. The FORMULA is a string containing combinations of
%   geometry object tags and operations ('+' union/join, '-' subtract,
%   '&' intersect). TOL (default 1e-5) specifies the tolerance for CSG
%   boundary operations. The status return code STAT is zero if the
%   operation has been applied sucessfully, or a positive integer
%   corresponding to the number of failed geometry object
%   combinations. The IS_WARN flag (default false) enables throwing
%   warnings if the geometry objects cannot be combined.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_apply_formula, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_apply_formula', varargin{:} );
if( ~nargout ), clear varargout; end
