function [ varargout ] = mesh_repair( varargin )
%MESH_REPAIR Repair 3D surface mesh.
%
%   [ P, F, N, IND, IER ] = MESH_REPAIR( P, F, N, IND, TOL_P ) Repair
%   3D surface mesh data including deduplication of nodes, and removal
%   of duplicated and collapsed facets. Input data is the vertex
%   coordinates P (n_p x 3), faces F (n_f x 3), and optional outward
%   directed normal vectors N (n_f x 3). IND is an optional facet
%   index vector to the original order.
%
%   If the optional vertex deduplication tolerance TOL_P is omitted,
%   the repair algorithm is run with successively increased tolerance
%   removing vertices until a water tight surface mesh is
%   achieved. IER indicates output status (0 for successful mesh
%   repair with a water tight mesh, and 1 otherwise).
%
%   Example:
%
%     offset = 0.1;
%     p = [ 0 0 0; 1 0 0; 0 0 1 ;
%           0 0 0; 0 1 0; 1 0 0 ;
%           0 0 0; 0 0 1; 0 1 0 ;
%           1+offset 0 0; 0 1+offset 0; 0 0 1+offset ];
%     f = [ 1:3; 4:6; 7:9; 10:12 ];
%
%     [p,f,n,ind,ier] = mesh_repair( p, f )
%
%   See also MESH_ANALYZE, DEDUPLICATE, POLYGON_NORMAL_AREA

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help mesh_repair, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'mesh_repair', varargin{:} );
if( ~nargout ), clear varargout; end
