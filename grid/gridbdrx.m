function [ varargout ] = gridbdrx( varargin )
%GRIDBDRX Reconstruct interior boundaries.
%
%   [ GRID ] = GRIDBDRX( GRID, ORDER ) Reconstructs interior
%   boundaries and updates the GRID.B field. ORDER (default true)
%   changes the order of boundaries so external boundaries are first,
%   and internal last.
%
%   See also GRIDBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridbdrx, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridbdrx', varargin{:} );
if( ~nargout ), clear varargout; end
