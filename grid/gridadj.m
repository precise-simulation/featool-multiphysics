function [ varargout ] = gridadj( varargin )
%GRIDADJ Generate grid cell adjacency information.
%
%   [ A ] = GRIDADJ( C, N_SDIM ) Generates an array A with indices to
%   neighboring cells and zeros where no neighbors are found. Input
%   data is an array C with cell connectivities and N_SDIM specifying
%   the number of space dimensions.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       c           (n_v,n_c)              Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers specifying the cells.
%       n_sdim      scalar                 Number of space dimensions.
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       a           (n_e,n_c)              Array with cell adjacency pointers where
%                                          the entries in each column correspond to
%                                          the neighbouring cell number for each
%                                          edge or face. A zero entry indicates that
%                                          the edge/face belongs to an external boundary.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridadj, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridadj', varargin{:} );
if( ~nargout ), clear varargout; end
