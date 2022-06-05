%SF_LINE_H3 Third order 1D C1 Hermite shape functions for lines.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_H3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates C1 Hermite shape functions on 1D line elements with value and
%   first derivatives defined in the nodes. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-4            Local basis function to evaluate
%       xi          array  [2,1]           Local coordinates of evaluation point
%       aInvJac     [n,3]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [2,4]                  Number of local degrees of freedom on
%                                          vertices, edges, faces, cell interiors,
%                                          and vertices without boundary conditions
%       xLDof       [2,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_LINE_P3

% Copyright 2013-2022 Precise Simulation, Ltd.


















