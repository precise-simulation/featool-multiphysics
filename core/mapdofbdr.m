function [ varargout ] = mapdofbdr( varargin )
%MAPDOFBDR Create degree of freedom and boundary maps.
%
%   [ N_GDOF, ADOFMAP, ABDRMAP ] = MAPDOFBDR( VARARGIN ) Computes a
%   connectivity map to the global degrees of freedom on each cell, and
%   optionally computes a boundary index vector for all degrees of freedom.
%   Input parameters can either be old type call (P,C,A,B,NLDOF) or
%   (C,B,NLDOF,XI) where the local dof coordinates XI is optional.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       p           (2,n_p)                Array with grid point coordinates
%       c           (n_vc,n_c)             Array with cell connectivities, points
%                                          to vertices in p for each cell
%       a           (n_cf,n_c)             Cell adjacency information, points to
%                                          adjacent cells for each edge. If the edge
%                                          is on a boundary the a entry is zero
%       b           (4-6,n_cb)             Boundary indicies for degress of freedom
%       nLDof       (n_bdgroups,4)         Number of local degrees of freedom on vertices
%                                          edges, and cell interiors for boundary groups
%       xi          (2-4,sum(nLDof))       Coordinates of local dofs (optional).
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       n_gdof      (scalar)               Total number of degrees of freedom
%       aDofMap     (n_ldof,n_c)           Connectivity map giving degree of freedom
%                                          numbers for each local dof on all cells
%       aBdrMap     (5+n_sdim,n_bdof)      Boundary dof map. First row is the cell number,
%                                          followed by edge/face, boundary, global and
%                                          local dof numbers, and local coordinates
%                                          on edges/faces.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help mapdofbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'mapdofbdr', varargin{:} );
if( ~nargout ), clear varargout; end
