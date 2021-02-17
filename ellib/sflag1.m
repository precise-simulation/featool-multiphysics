function [ vBase, nLDof, xLDof, sfun ] = sflag1( i_eval, n_sdim, n_vert, varargin )
%SFLAG1 Shape function driver (first order P1/Q1 Lagrange polynomials).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SFLAG1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming linear P1 Lagrange shape functions with values
%   defined in the nodes.
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
%   See also SF_SIMP_P1, SF_QUAD_Q3, SF_HEX_Q3

% Copyright 2013-2021 Precise Simulation, Ltd.


if ( (n_vert==n_sdim+1) || (n_sdim==1) )   % Simplex cells (1D line, 2D triangle, 3D tetrahedron).

  [vBase,nLDof,xLDof,sfun] = sf_simp_P1(i_eval,n_sdim,n_vert,varargin{:});

else   % Quadrilateral and hexahedral cells.

  switch n_sdim

    case 1

      error('sflag1: 1D only supports barycentric coordinates.')

    case 2   % 2D Quadrilateral cells.

      [vBase,nLDof,xLDof,sfun] = sf_quad_Q1(i_eval,n_sdim,n_vert,varargin{:});

    case 3   % 3D Hedahedral cells.

      [vBase,nLDof,xLDof,sfun] = sf_hex_Q1(i_eval,n_sdim,n_vert,varargin{:});

  end
end
