%SFLAG5 Shape function driver (fifth order P5/Q5 Lagrange polynomials).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SFLAG5( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming quintic Lagrange shape functions with values
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
%   See also SF_LINE_P5, SF_TRI_P5, SF_TET_P5, SF_QUAD_Q5, SF_HEX_Q5

% Copyright 2013-2025 Precise Simulation, Ltd.


















