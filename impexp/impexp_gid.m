function [ varargout ] = impexp_gid( varargin )
%IMPEXP_GID Import/export grid in GID msh format.
%
%   [ GRID ] = IMPEXP_GID( FILE_NAME, MODE, GRID, FID_LOG ) Import or export of
%   a GiD mesh ASCII format .msh file. FILE_NAME is a string specifying the (root)
%   file name to process. MODE can either be a string indicating import (no boundary
%   reconstruction), import_bdr (boundary reconstruction with gridbdr), or export.
%   GRID is a grid struct input for export and output for import. FID_LOG is an
%   optional log file handle for message output (negative for gui output or empty
%   for no output).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_gid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_gid', varargin{:} );
if( ~nargout ), clear varargout; end
