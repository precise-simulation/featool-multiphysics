function [ varargout ] = geom_revolve_face( varargin )
%GEOM_REVOLVE_FACE Revolve face to solid.
%
%   [ SOUT, NEW_TAG ] = GEOM_REVOLVE_FACE( SIN, TAG, FACE, TH, P, V, PIT )
%   Revolve FACE from geometry object with TAG an angle TH
%   (degrees) around axis defined by point P and direction vector V.
%   PIT is optionally the pitch (distance per revolution) to offset
%   the revolution in the axial direction (to form a coil/helix).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_revolve_face, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_revolve_face', varargin{:} );
if( ~nargout ), clear varargout; end
