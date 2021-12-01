function [ varargout ] = gridedge( varargin )
%GRIDEDGE Generate grid cell edge numbering.
%
%   [ E, EV ] = GRIDEDGE( C, N_SDIM ) Generates an array E with global
%   edge numbers. Input data is an array C with cell connectivities
%   and N_SDIM specifying the number of space dimensions.
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
%       e           (n_e,n_c)              Array with global edge numbers for all cells.
%       ev          (n_e,2)                Array with indices to vertices for each edge.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridedge, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridedge', varargin{:} );
if( ~nargout ), clear varargout; end
