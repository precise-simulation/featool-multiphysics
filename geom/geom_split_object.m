function [ varargout ] = geom_split_object( varargin )
%GEOM_SPLIT_OBJECT Split geometry object(s) by plane or line.
%
%   [ SOUT, NEW_TAGS ] = GEOM_SPLIT_OBJECT( SIN, TAGS, P, V ) Splits
%   the geometry object(s) with TAGS by the plane defined by the point P
%   and normal vector V in 3D (or line in 2D with V being the
%   direction vector). SIN/SOUT is either a geometry or fea struct.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_split_object, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_split_object', varargin{:} );
if( ~nargout ), clear varargout; end
