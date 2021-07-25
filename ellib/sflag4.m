function [ vBase, nLDof, xLDof, sfun ] = sflag4( i_eval, n_sdim, n_vert, varargin )
%SFLAG4 Shape function driver (fourth order P4/Q4 Lagrange polynomials).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SFLAG4( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming quartic Lagrange shape functions with values
%   defined in the nodes, and edges, (and also faces, and cell centers for
%   quadrilaterals, tetrahedra, and hexahedrals).
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1-3            Number of space dimensions
%       n_vert      scalar: 2-8            Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim(+1)]           Local coordinates of evaluation point
%       aInvJac     [n,n_sdim(+1)*n_sdim]  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [n_sdim,n_ldof]        Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_LINE_P4, SF_TRI_P4, SF_TET_P4, SF_QUAD_Q4, SF_HEX_Q4

% Copyright 2013-2021 Precise Simulation, Ltd.


if ( (n_vert==n_sdim+1) || (n_sdim==1) )   % Simplex cells (1D line, 2D triangle, 3D tetrahedron).

  switch n_sdim

    case 1   % 1D Line cells.

      [vBase,nLDof,xLDof,sfun] = sf_line_P4(i_eval,n_sdim,n_vert,varargin{:});

    case 2   % 2D Triangular cells.

      [vBase,nLDof,xLDof,sfun] = sf_tri_P4(i_eval,n_sdim,n_vert,varargin{:});

    case 3   % 3D Tetrahedral cells.

      [vBase,nLDof,xLDof,sfun] = sf_tet_P4(i_eval,n_sdim,n_vert,varargin{:});

  end

else

  switch n_sdim

    case 2   % 2D Quadrilateral cells.

      [vBase,nLDof,xLDof,sfun] = sf_quad_Q4(i_eval,n_sdim,n_vert,varargin{:});

    case 3   % 3D Hedahedral cells.

      [vBase,nLDof,xLDof,sfun] = sf_hex_Q4(i_eval,n_sdim,n_vert,varargin{:});

  end

end
