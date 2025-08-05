%SF_SIMP_RT1 Linear vector (Raviart-Thomas) divergence shape function for simplices.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_RT1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates linear vector Raviart-Thomas divergence shape functions on simplices with
%   values defined in the nodes. XI is Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 2/3            Number of space dimensions
%       n_vert      scalar: 3/4            Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim+1]             Local coordinates of evaluation point
%       aInvJac     [n,n_sdim+1*n_sdim]    Inverse of transformation Jacobian
%       vBase       [n,1,2/3]              Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n,1,2/3]              Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [n_sdim,n_ldof]        Local coordinates of local dofs
%       sfun        string                 Function name of called shape function

% Copyright 2013-2025 Precise Simulation, Ltd.

























