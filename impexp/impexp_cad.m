function [ varargout ] = impexp_cad( varargin )
%IMPEXP_CAD Import/export geometry in CAD formats.
%
%   [ GEOM ] = IMPEXP_STEP( FILE_NAME, MODE, GEOM, FORMAT ) Import and
%   export of 3D geometries and in BREP, and IGES, STEP CAD formats.
%
%   FILE_NAME is a string which specifies the file name to
%   process.
%
%   MODE is a string indicating either IMPORT (default) or EXPORT.
%
%   During export, a FEATool geometry struct GEOM must also be
%   provided. If provided during import the output geometry will
%   include the original geometry objects.
%
%   FORMAT is a string designating CAD format to import/export (valid
%   formats are brep, iges, and step, default taken from the input
%   file extension).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_cad, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_cad', varargin{:} );
if( ~nargout ), clear varargout; end
