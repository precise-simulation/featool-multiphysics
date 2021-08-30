function [ varargout ] = impexp_vtk( varargin )
%EXPORT_VTK Export grid and data in VTK format.
%
%   [ GRID ] = IMPEXP_VTK( FILENAME, MODE, DATA, EXPR, SOLNUM )
%
%   Import of grid data and export of grid points and corresponding
%   solution/expression data in VTK 3.0 ASCII format (non XML).
%
%   FILENAME specifies the input/output file name. MODE can either be
%   a string indicating import (no boundary reconstruction),
%   import_bdr (boundary reconstruction with gridbdr), or export
%   (default). DATA must be a valid grid or FEA data struct. EXPR and
%   SOLNUM are optional arguments for export, where EXPR is a cell
%   array with expressions to process. Expressions will be evaluated
%   in point corresponding to SCALARS/POINT_DATA format, and dependent
%   variables present in the FEA struct will be prepended to
%   EXPR. SOLNUM is optionally the solution number of which to
%   evaluate the expressions. FID_LOG is an optional log file handle
%   for message output (negative for gui output or empty for no
%   output).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_vtk, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_vtk', varargin{:} );
if( ~nargout ), clear varargout; end
