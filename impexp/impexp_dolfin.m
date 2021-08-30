function [ varargout ] = impexp_dolfin( varargin )
%IMPEXP_DOLFIN Import/export grid in FEniCS/Dolfin XML format.
%
%   [ GRID ] = IMPEXP_DOLFIN( FILE_NAME, MODE, DATA, USE_MESHFCN, FID )
%   Import or export of a FEniCS/Dolfin (ASCII .xml) grid and data
%   format. FILE_NAME is a string specifying the (root) file name to
%   process. MODE can either be a string indicating import (no
%   boundary reconstruction), import_bdr (boundary reconstruction with
%   gridbdr), or export. For export, DATA can be either a full fea
%   struct or just the grid struct. USE_MESHFCN is a boolean flag to
%   write subdomain numbers in a mesh function instead of the
%   (default) domain field. A GRID struct is output when importing
%   FID is an optional log file handle for message output
%   (negative for gui output or empty for no output).
%
%   Example:
%
%      1) Export and re-import of a unit square with a hole grid.
%
%      grid1 = quad2tri(holegrid());
%      impexp_hdf5( 'featool-fenics-mesh.xml', 'export', grid1 )
%      grid2 = impexp_hdf5( 'featool-fenics-mesh.h5', 'import' )
%      subplot(1,2,1), plotgrid(grid1), title('grid1')
%      subplot(1,2,2), plotgrid(grid2), title('grid2')
%      is_ok = gridcheck(grid2) == 0
%
%   See also FENICS, FENICS_IMPORT, IMPEXP_HDF5

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_dolfin, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_dolfin', varargin{:} );
if( ~nargout ), clear varargout; end
