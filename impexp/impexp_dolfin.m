function [ varargout ] = impexp_dolfin( varargin )
%IMPEXP_DOLFIN Import/export grid in Dolfin XML format.
%
%   [ GRID ] = IMPEXP_DOLFIN( FILE_NAME, MODE, DATA, USE_MESHFCN, FID_LOG ) Import
%   or export of a Dolfin/fenics (.xml) grid and data format. FILE_NAME is a string
%   specifying the (root) file name to process. MODE can either be a string indicating
%   import (no boundary reconstruction), import_bdr (boundary reconstruction with gridbdr),
%   or export. For export, DATA can be either a full fea struct or just the grid
%   struct. USE_MESHFCN is a boolean flag to write subdomain numbers in a mesh function
%   instead of the (default) domain field. A GRID struct is output when importing.
%   FID_LOG is an optional log file handle for message output (negative for gui output
%   or empty for no output).

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_dolfin, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_dolfin', varargin{:} );
if( ~nargout ), clear varargout; end
