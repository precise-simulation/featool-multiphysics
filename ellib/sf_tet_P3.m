%SF_TET_P3 Third order Lagrange shape functions for tetrahedra (P3).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_TET_P3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming third order P3 Lagrange shape functions on 3D tetrahedral elements
%   with values defined in the nodes, edges, and faces. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 3              Number of space dimensions
%       n_vert      scalar: 4              Number of vertices per cell
%       i_dof       scalar: 1-16           Local basis function to evaluate
%       xi          array [4,1]            Local coordinates of evaluation point
%       aInvJac     [n,12]                 Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [4,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SFLAG3

% Copyright 2013-2024 Precise Simulation, Ltd.












































