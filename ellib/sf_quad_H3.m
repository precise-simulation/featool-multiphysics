%SF_QUAD_H3 Third order 2D C1 Hermite shape function for quadrilaterals.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_QUAD_H3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates C1 Hermite shape functions on quadrilaterals with values defined in the nodes,
%   cell center. XI are [-1..1]^2 reference coordinates. (Limited to rectangular grid cells.)
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  2             Number of space dimensions
%       n_vert      scalar:  4             Number of vertices per cell
%       i_dof       scalar: 1-16           Local basis function to evaluate
%       xi          [n_sdim]               Local coordinates of evaluation point
%       aInvJac     [n,n_sdim*n_sdim]      Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [3,4]                  Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [n_sdim,n_ldof]        Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_QUAD_Q3

% Copyright 2013-2023 Precise Simulation, Ltd.



















