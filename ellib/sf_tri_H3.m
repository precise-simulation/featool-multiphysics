function [ vBase, nLDof, xLDof, sfun ] = sf_tri_H3( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
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

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [3 0 0 1;
         3 0 0 0;
         3 0 0 0];
xLDof = [repmat([1 0 0;
                 0 1 0;
                 0 0 1],1,3) [1/3;1/3;1/3]];
sfun  = 'sf_tri_H3';


switch i_eval

  case 1

    switch i_dof
      case 1
        vBase = - 2*xi(1)^3 + 3*xi(1)^2 - 7*xi(2)*xi(3)*xi(1);
      case 2
        vBase = - 2*xi(2)^3 + 3*xi(2)^2 - 7*xi(1)*xi(3)*xi(2);
      case 3
        vBase = - 2*xi(3)^3 + 3*xi(3)^2 - 7*xi(1)*xi(2)*xi(3);
      case 4
        vBase = - xi(1)*xi(2)*((-aInvJac(:,6).*aInvJac(:,7)))*(2*xi(1) + xi(2) - 1) - xi(1)*xi(3)*(aInvJac(:,5).*aInvJac(:,7))*(2*xi(1) + xi(3) - 1);
      case 5
        vBase = xi(1)*xi(2)*((-aInvJac(:,6).*aInvJac(:,7)))*(xi(1) + 2*xi(2) - 1) - xi(2)*xi(3)*((-aInvJac(:,4).*aInvJac(:,7)))*(2*xi(2) + xi(3) - 1);
      case 6
        vBase = xi(1)*xi(3)*(aInvJac(:,5).*aInvJac(:,7))*(xi(1) + 2*xi(3) - 1) + xi(2)*xi(3)*((-aInvJac(:,4).*aInvJac(:,7)))*(xi(2) + 2*xi(3) - 1);
      case 7
        vBase = - xi(1)*xi(2)*(aInvJac(:,3).*aInvJac(:,7))*(2*xi(1) + xi(2) - 1) - xi(1)*xi(3)*((-aInvJac(:,2).*aInvJac(:,7)))*(2*xi(1) + xi(3) - 1);
      case 8
        vBase = xi(1)*xi(2)*(aInvJac(:,3).*aInvJac(:,7))*(xi(1) + 2*xi(2) - 1) - xi(2)*xi(3)*(aInvJac(:,1).*aInvJac(:,7))*(2*xi(2) + xi(3) - 1);
      case 9
        vBase = xi(1)*xi(3)*((-aInvJac(:,2).*aInvJac(:,7)))*(xi(1) + 2*xi(3) - 1) + xi(2)*xi(3)*(aInvJac(:,1).*aInvJac(:,7))*(xi(2) + 2*xi(3) - 1);
      case 10
        vBase = 27*xi(1)*xi(2)*xi(3);
    end

  case {2,3}

    switch i_dof

      case 1
        dNdxi1 = - 6*xi(1)^2 + 6*xi(1) - 7*xi(2)*xi(3);
        dNdxi2 = -7*xi(1)*xi(3);
        dNdxi3 = -7*xi(1)*xi(2);
      case 2
        dNdxi1 = -7*xi(2)*xi(3);
        dNdxi2 = - 6*xi(2)^2 + 6*xi(2) - 7*xi(1)*xi(3);
        dNdxi3 = -7*xi(1)*xi(2);
      case 3
        dNdxi1 = -7*xi(2)*xi(3);
        dNdxi2 = -7*xi(1)*xi(3);
        dNdxi3 = - 6*xi(3)^2 + 6*xi(3) - 7*xi(1)*xi(2);
      case 4
        dNdxi1 = - xi(2)*((-aInvJac(:,6).*aInvJac(:,7)))*(2*xi(1) + xi(2) - 1) - xi(3)*(aInvJac(:,5).*aInvJac(:,7))*(2*xi(1) + xi(3) - 1) - 2*xi(1)*xi(2)*((-aInvJac(:,6).*aInvJac(:,7))) - 2*xi(1)*xi(3)*(aInvJac(:,5).*aInvJac(:,7));
        dNdxi2 = -xi(1)*((-aInvJac(:,6).*aInvJac(:,7)))*(2*xi(1) + 2*xi(2) - 1);
        dNdxi3 = -xi(1)*(aInvJac(:,5).*aInvJac(:,7))*(2*xi(1) + 2*xi(3) - 1);
      case 5
        dNdxi1 = xi(2)*((-aInvJac(:,6).*aInvJac(:,7)))*(2*xi(1) + 2*xi(2) - 1);
        dNdxi2 = xi(1)*((-aInvJac(:,6).*aInvJac(:,7)))*(xi(1) + 2*xi(2) - 1) - xi(3)*((-aInvJac(:,4).*aInvJac(:,7)))*(2*xi(2) + xi(3) - 1) + 2*xi(1)*xi(2)*((-aInvJac(:,6).*aInvJac(:,7))) - 2*xi(2)*xi(3)*((-aInvJac(:,4).*aInvJac(:,7)));
        dNdxi3 = -xi(2)*((-aInvJac(:,4).*aInvJac(:,7)))*(2*xi(2) + 2*xi(3) - 1);
      case 6
        dNdxi1 = xi(3)*(aInvJac(:,5).*aInvJac(:,7))*(2*xi(1) + 2*xi(3) - 1);
        dNdxi2 = xi(3)*((-aInvJac(:,4).*aInvJac(:,7)))*(2*xi(2) + 2*xi(3) - 1);
        dNdxi3 = xi(1)*(aInvJac(:,5).*aInvJac(:,7))*(xi(1) + 2*xi(3) - 1) + xi(2)*((-aInvJac(:,4).*aInvJac(:,7)))*(xi(2) + 2*xi(3) - 1) + 2*xi(1)*xi(3)*(aInvJac(:,5).*aInvJac(:,7)) + 2*xi(2)*xi(3)*((-aInvJac(:,4).*aInvJac(:,7)));
      case 7
        dNdxi1 = - 2*xi(1)*xi(2)*(aInvJac(:,3).*aInvJac(:,7)) - 2*xi(1)*xi(3)*((-aInvJac(:,2).*aInvJac(:,7))) - xi(2)*(aInvJac(:,3).*aInvJac(:,7))*(2*xi(1) + xi(2) - 1) - xi(3)*((-aInvJac(:,2).*aInvJac(:,7)))*(2*xi(1) + xi(3) - 1);
        dNdxi2 = -xi(1)*(aInvJac(:,3).*aInvJac(:,7))*(2*xi(1) + 2*xi(2) - 1);
        dNdxi3 = -xi(1)*((-aInvJac(:,2).*aInvJac(:,7)))*(2*xi(1) + 2*xi(3) - 1);
      case 8
        dNdxi1 = xi(2)*(aInvJac(:,3).*aInvJac(:,7))*(2*xi(1) + 2*xi(2) - 1);
        dNdxi2 = 2*xi(1)*xi(2)*(aInvJac(:,3).*aInvJac(:,7)) - 2*xi(2)*xi(3)*(aInvJac(:,1).*aInvJac(:,7)) + xi(1)*(aInvJac(:,3).*aInvJac(:,7))*(xi(1) + 2*xi(2) - 1) - xi(3)*(aInvJac(:,1).*aInvJac(:,7))*(2*xi(2) + xi(3) - 1);
        dNdxi3 = -xi(2)*(aInvJac(:,1).*aInvJac(:,7))*(2*xi(2) + 2*xi(3) - 1);
      case 9
        dNdxi1 = xi(3)*((-aInvJac(:,2).*aInvJac(:,7)))*(2*xi(1) + 2*xi(3) - 1);
        dNdxi2 = xi(3)*(aInvJac(:,1).*aInvJac(:,7))*(2*xi(2) + 2*xi(3) - 1);
        dNdxi3 = 2*xi(1)*xi(3)*((-aInvJac(:,2).*aInvJac(:,7))) + 2*xi(2)*xi(3)*(aInvJac(:,1).*aInvJac(:,7)) + xi(1)*((-aInvJac(:,2).*aInvJac(:,7)))*(xi(1) + 2*xi(3) - 1) + xi(2)*(aInvJac(:,1).*aInvJac(:,7))*(xi(2) + 2*xi(3) - 1);
      case 10
        dNdxi1 = 27*xi(2)*xi(3);
        dNdxi2 = 27*xi(1)*xi(3);
        dNdxi3 = 27*xi(1)*xi(2);
    end

    if( i_eval==2 )
      vBase = aInvJac(:,1).*dNdxi1 + aInvJac(:,2).*dNdxi2 + aInvJac(:,3).*dNdxi3;
    else
      vBase = aInvJac(:,4).*dNdxi1 + aInvJac(:,5).*dNdxi2 + aInvJac(:,6).*dNdxi3;
    end

  otherwise
    vBase = 0;

end
