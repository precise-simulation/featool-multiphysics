function [ vBase, nLDof, xLDof, sfun ] = sf_simp_P1( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_SIMP_P1 Linear Lagrange shape function for simplices (P1).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_P1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming linear P1 Lagrange shape functions on simplices with
%   values defined in the nodes. XI is Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1-3            Number of space dimensions
%       n_vert      scalar: 2-4            Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim+1]             Local coordinates of evaluation point
%       aInvJac     [n,n_sdim+1*n_sdim]    Inverse of transformation Jacobian
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
%   See also SFLAG1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [n_vert 0 0 0];
xLDof = eye(n_vert);
sfun  = 'sf_simp_P1';


switch i_eval    % Evaluation type flag.

  case 1         % Evaluation of function values.

    vBase = xi(i_dof);

  case {2,3,4}   % Evaluation of first derivatives.

    i_col = i_dof+(i_eval-2)*n_vert;

    vBase = aInvJac(:,i_col);

  otherwise
    vBase = 0;

end
