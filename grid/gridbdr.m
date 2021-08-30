function [ varargout ] = gridbdr( varargin )
%GRIDBDR Generate grid boundary information.
%
%   [ B ] = GRIDBDR( P, C, A, THLIM, FIX_BDR, B0 ) Given arrays P, C and A with
%   cell connectivities and adjacency information generates an array with
%   boundary information. The last three input arguments are optional.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       p           (n_sdim,n_p)           Coordinate array for grid vertices.
%       c           (n_pcell,n_cells)      Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers.
%       a           (n_pcell,n_cells)      Array with cell adjacency pointers where
%                                          the entries in each column correspond to
%                                          the neighbouring cell number for each edge
%                                          or face. A zero entry indicates that the
%                                          edge/face belongs to a boundary.
%       thlim       scalar {0.77}          boundary edge detection limit, splits boundaries
%                                          with cos(normal angles) less than thlim.
%       fix_bdr     logical {0}            collect all boundary edges to the same boundary
%       b0          (4/5/6,n_b0)           Array of partially precomputed boundaries
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       b           (4/5/6,n_b)            Array with boundary information, the first row
%                                          indicates the cell number, followed by cell
%                                          edge/face, boundary number, and the components
%                                          of the normal vector for each boundary edge/face
%
%   See also GRIDBDRE, GRIDBDRX

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridbdr', varargin{:} );
if( ~nargout ), clear varargout; end
