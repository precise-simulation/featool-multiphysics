function [ varargout ] = findbdr( varargin )
%FINDBDR Find boundary indices for an expression.
%
%   [ BDR, INDB ] = FINDBDR( SIN, S_EXPR, MATCH_ALL ) Evaluates S_EXPR
%   for boundary nodes in the grid or fea struct and returns the boundary
%   numbers BDR (corresponding to values in grid.b(3,:)) for which the
%   expression is true, and also indices INDB to the corresponding boundary
%   segments (columns in grid.b). The flag MATCH_ALL requires all nodes on
%   a boundary to fulfill S_EXPR to return the found indices (default). If
%   MATCH_ALL is false all indices to all boundaries BDR, edges, and faces
%   INDB are returned.
%
%   Examples:
%
%      1) Set boundary number 2 on all boundary edges/faces where x>0.5.
%
%      [~,indb]  = findbdr( grid, 'x>0.5', false );
%      b(3,indb) = 2;

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help findbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'findbdr', varargin{:} );
if( ~nargout ), clear varargout; end
