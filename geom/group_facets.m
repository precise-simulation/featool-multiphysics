function [ varargout ] = group_facets( varargin )
%GROUP_FACETS Grouping of triangular 3D facets.
%
%   [ B ] = GROUP_FACETS( P, F, N, TOL_N, N_F_GRP, MERGE_SINGLE, THLIM )
%   Collect 3D triangles according to groups of connected facets with
%   aligned normals and returns a geometry object boundary struct B.
%
%   TOL_N controls the tolerance for grouping of facets
%   into boundaries by comparing the difference between the normals
%   (default 1e-2 relative tolerance). Setting this to a large value
%   effectively disables facet grouping.
%
%   N_F_GRP is an integer specifying the maximum number of facets in
%   groups for which to merge together to new boundaries (default 1
%   only merging ungrouped facets).
%
%   MERGE_SINGLE toggles merging of single ungrouped facets (default true).
%
%   THLIM Sorts facets (with n_facets <= N_F_GRP) in to groups with
%   normals within THLIM degrees (default 0 to skip grouping).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help group_facets, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'group_facets', varargin{:} );
if( ~nargout ), clear varargout; end
