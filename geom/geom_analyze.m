function [ varargout ] = geom_analyze( varargin )
%GEOM_ANALYZE Analyze geometry objects.
%
%   [ SOUT ] = GEOM_ANALYZE( SIN, TOL, IS_WARN, USE_GEOMTOOL, IS_SPLIT )
%   Analyzes geometry (objects) SIN and returns a decomposed geometry
%   with non-overlapping geometry objects. TOL (default 1e-5)
%   specifies the tolerance for CSG boundary operations. The IS_WARN
%   flag (default false) enables throwing warnings if the geometry
%   objects cannot be combined. The USE_GEOMTOOL flag (default true)
%   enables using the GEOMTool engine for geometry decomposition if
%   available. SPLIT (default true) returns separate objects, or a
%   single combined one (false, only with geomtool).
%
%   See also GEOMCFG

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_analyze, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_analyze', varargin{:} );
if( ~nargout ), clear varargout; end
