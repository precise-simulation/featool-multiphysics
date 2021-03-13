function [ varargout ] = gridface( varargin )
%GRIDFACE Generate grid cell face numbering.
%
%   [ F, FV ] = GRIDFACE( C, B ) Generates an array F with global face
%   numbers and FV with grid point indices per face. Input data is an
%   array C with cell connectivities, and optionally the boundary
%   struct B if upper triangular (OpenFOAM) ordering should be used
%   (default normal ordering).
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       c           (n_v,n_c)              Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers specifying the cells.
%       b           struct                 Boundary struct to use upper triangular ordering
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (n_f,n_c)              Array with global face numbers for all cells.
%       fv          (n_f,3/4)              Array with indices to vertices for each face.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridface, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridface', varargin{:} );
if( ~nargout ), clear varargout; end
