%SFNONCON1 Shape function driver (first order P1/Q1 Lagrange polynomials).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SFNONCON1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates nonconforming linear (P-1/Q1~) shape functions with degrees
%   of freedom defined on the cell edges/faces.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 2-3            Number of space dimensions
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
%   See also SF_SIMP_P1NC, SF_QUAD_Q1NC, SF_HEX_Q1NC

% Copyright 2013-2025 Precise Simulation, Ltd.












