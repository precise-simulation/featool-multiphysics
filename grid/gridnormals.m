function [ varargout ] = gridnormals( varargin )
%GRIDNORMALS Compute normal vectors.
%
%   [ N ] = GRIDNORMALS( P, C, B ) Computes an array of normal
%   vectors. P is an array of grid point coordinates and C grid
%   vertex connectivities. B specify boundary edges/faces with
%   the fist row B(1,:) indicating cell numbers and the second
%   B(2,:) local cell edge/face numbering.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridnormals, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridnormals', varargin{:} );
if( ~nargout ), clear varargout; end
