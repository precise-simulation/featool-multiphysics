%SF_TRI_H3 Third order 2D C1 Hermite shape functions for triangles.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_TRI_H3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates C1 Hermite shape functions on quadrilaterals with values defined in the nodes,
%   and cell center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 2              Number of space dimensions
%       n_vert      scalar: 3              Number of vertices per cell
%       i_dof       scalar: 1-10           Local basis function to evaluate
%       xi          array [3,1]            Local coordinates of evaluation point
%       aInvJac     [n,6]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [3,4]                  Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [3,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_TRI_P1

% Copyright 2013-2025 Precise Simulation, Ltd.












