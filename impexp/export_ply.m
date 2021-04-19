function [ varargout ] = export_ply( varargin )
%EXPORT_PLY Export geometry in PLY format.
%
%   EXPORT_PLY( FILE_NAME, GEOM ) Export of 3D CAD geometries in ASCII
%   Polygon File Format/Stanford Triangle Format. FILE_NAME is a
%   string which specifies the file name to process. GEOM is the
%   geometry to export.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help export_ply, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'export_ply', varargin{:} );
if( ~nargout ), clear varargout; end
