function [ varargout ] = gridbdre( varargin )
%GRIDBDRE Generate grid edge boundary information.
%
%   [ BE, E, EV ] = GRIDBDRE( B, C ) Given arrays B and C with cell
%   connectivities and face boundary information generates an array
%   with edge boundary information.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       b           (3/6,n_b)              Array with boundary information, the first row
%                                          indicates the cell number, followed by cell
%                                          face, boundary number, and optionally the
%                                          normal vectors for each boundary face
%       c           (n_pcell,n_cells)      Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers.
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       be          (5,n_b)                Array with boundary information, the first row
%                                          indicates the global cell number, followed by
%                                          local cell face, local edge, and global edge
%                                          and edge boundary numbers.
%       e           (n_e,n_c)              Array with global edge numbers for all cells.
%       ev          (n_e,2)                Array with indices to vertices for each edge.
%
%   See also GRIDBDR, GRIDEDGE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridbdre, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridbdre', varargin{:} );
if( ~nargout ), clear varargout; end
