function [ varargout ] = impexp_gmv( varargin )
%IMPEXP_GMV Import/export grid in GMV format.
%
%   [ GRID ] = IMPEXP_GMV( FILE_NAME, MODE, DATA, N_SDIM, FID_LOG, SOLNUM )
%   Import or export of a General Mesh Viewer (GMV) .gmv format data file.
%
%   FILE_NAME is a string specifying the (root) file name to
%   process. MODE can either be a string indicating import (no
%   boundary reconstruction), import_bdr (boundary reconstruction with
%   gridbdr), or export. For export, DATA can be either a whole data
%   struct or just the grid struct. N_SDIM is the number of space
%   dimensions. And SOLNUM is an optional solution number to export. A
%   GRID struct is the returned output after importing. FID_LOG is an
%   optional log file handle for message output (negative for gui
%   output or empty for no output). For export, expressions stored in
%   a DATA.post.expr cell array will be evaluated and exported as scalars.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_gmv, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_gmv', varargin{:} );
if( ~nargout ), clear varargout; end
