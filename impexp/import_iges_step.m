function [ varargout ] = import_iges_step( varargin )
%IMPEXP_IGES_STEP Import IGES/STEP CAD file.
%
%   [ GEOM ] = IMPORT_IGES_STEP( FILE_NAME ) Import IGES or STEP CAD
%   file and converts to a geom struct. FILE_NAME is a string
%   specifying the IGES/STEP file to read.

% Copyright 2013-2019 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help import_iges_step, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'import_iges_step', varargin{:} );
if( ~nargout ), clear varargout; end
