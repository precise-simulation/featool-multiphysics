function [ varargout ] = gridvert( varargin )
%GRIDVERT Generate grid vertex adjacency map.
%
%   [ V, VV ] = GRIDVERT( C, INDC, N_SDIM ) Generates an sparse map V
%   which indicates which cells are connected to each other through
%   the nodes.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       c           (n_v,n_c)              Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers specifying the cells.
%       indc        vector                 Index to cells for which to compute map
%                                          (all cells by default).
%       n_sdim      scalar                 Optional input argument indicating computation.
%                                          of connectivities for connected edges.
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       v           sparse(n_c,n_c)        Sparse array indicating which grid cells are
%                                          connected through the grid points. Or with
%                                          three input arguments a which grid cells are
%                                          connected though the edges (incl. n_sdim).
%       vv          sparse(n_p,n_p)        Sparse array which for each row indicates
%                                          which nodes/edges are connected through cells.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridvert, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridvert', varargin{:} );
if( ~nargout ), clear varargout; end
