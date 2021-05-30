function [ varargout ] = impexp_hdf5( varargin )
%IMPEXP_HDF5 Import/export grid FEniCS HDF5 format.
%
%   [ GRID ] = IMPEXP_HDF5( FILE_NAME, MODE, DATA, FID ) Import or
%   export a FEniCS grid in hdf5 format. FILE_NAME is a string
%   specifying the file to process. MODE can be a string indicating
%   IMPORT (using automatic boundary reconstruction with gridbdr)
%   which returns a FEATool grid struct. If mode is EXPORT, the grid or
%   corresponding fea problem struct to export must be given in DATA.
%   FID is an optional log file handle for message output (negative
%   for gui output or empty for no output).
%
%   The following hdf5 datasets/fields are read/written:
%
%       dataset                  grid entry       size/type
%       -------------------------------------------------------------
%       /mesh/cell_indices       0:n_c-1          n_c/int64
%       /mesh/coordinates        grid.p           [n_sdim,n_p]/double
%       /mesh/topology           sort(grid.c)-1   [n_v,n_c]/int64
%       /mesh/values             grid.s-1         n_c/uint64
%
%       /boundary/coordinates    grid.p           [n_sdim,n_p]/double
%       /boundary/topology       gridedge/face    [n_ef,n_c]/int64
%       /boundary/values         grid.b(3,:)      n_b/uint64
%
%       /edge/coordinates        grid.p           [3,n_p]/double
%       /edge/topology           gridbdre         [2,n_c]/int64
%       /edge/values             be(5,:)          n_e/uint64
%       (only exported if DATA contains an .edg field)
%
%   Example:
%
%      1) Export and re-import of a unit circle grid.
%
%      grid1 = quad2tri(circgrid());
%      impexp_hdf5( 'featool-fenics-mesh.h5', 'export', grid1 )
%      grid2 = impexp_hdf5( 'featool-fenics-mesh.h5', 'import' )
%      subplot(1,2,1), plotgrid(grid1), title('grid1')
%      subplot(1,2,2), plotgrid(grid2), title('grid2')
%      is_ok = gridcheck(grid2) == 0
%
%   See also FENICS, FENICS_IMPORT, IMPEXP_DOLFIN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_hdf5, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_hdf5', varargin{:} );
if( ~nargout ), clear varargout; end
