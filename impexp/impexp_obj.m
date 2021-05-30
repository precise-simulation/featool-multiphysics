function [ varargout ] = impexp_obj( varargin )
%IMPEXP_OBJ Import/export geometry in OBJ format.
%
%   [ GEOM ] = IMPEXP_OBJ( FILE_NAME, MODE, GEOM, TEST2D ) Import and
%   export of 3D CAD geometries in Wavefront OBJ format. Only
%   triangular facets with the Obj fields v for vertices, f for faces,
%   and g for boundary groups are supported. The TEST2D flag (default
%   true) checks if all abs(z coordinates)<eps and reconstructs a 2D
%   planar geometry object.
%
%   FILE_NAME is a string which specifies the file name to
%   process.
%
%   MODE is a string indicating either IMPORT (default) or EXPORT.
%
%   During export, a FEATool geometry struct GEOM must also be
%   provided. If provided during import the output geometry will
%   include the original geometry objects.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_obj, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_obj', varargin{:} );
if( ~nargout ), clear varargout; end
