%SF_LINE_P5 1D Fifth order Lagrange shape functions for lines (P5).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_P5( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fifth order P5 Lagrange shape functions on 1D line elements
%   with values defined in the nodes and center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-6            Local basis function to evaluate
%       xi          array [2,1]            Local coordinates of evaluation point
%       aInvJac     [n,3]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [2,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_LINE_P1

% Copyright 2013-2025 Precise Simulation, Ltd.

















