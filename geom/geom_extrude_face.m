function [ varargout ] = geom_extrude_face( varargin )
%GEOM_EXTRUDE_FACE Extrude face to solid.
%
%   [ SOUT, NEW_TAG ] = GEOM_EXTRUDE_FACE( SIN, TAG, FACE, D, V )
%   Extrude FACE from geometry object with TAG a distance D in the
%   direction of the vector V.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_extrude_face, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_extrude_face', varargin{:} );
if( ~nargout ), clear varargout; end
