function [ vBase, nLDof, xLDof, sfun ] = sf_tet_P3( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
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

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [4 12 4 0];
xLDof = [1 0 0 0  2/3   0 1/3 2/3   0   0  1/3   0 2/3 1/3   0   0  1/3 1/3   0 1/3;
         0 1 0 0  1/3 2/3   0   0 2/3   0  2/3 1/3   0   0 1/3   0  1/3 1/3 1/3   0;
         0 0 1 0    0 1/3 2/3   0   0 2/3    0 2/3 1/3   0   0 1/3  1/3   0 1/3 1/3;
         0 0 0 1    0   0   0 1/3 1/3 1/3    0   0   0 2/3 2/3 2/3    0 1/3 1/3 1/3];
sfun  = 'sf_tet_P3';


switch i_eval

  case 1
    switch i_dof
      case 1
        vBase = (xi(1)*(3*xi(1) - 1)*(3*xi(1) - 2))/2;
      case 2
        vBase = (xi(2)*(3*xi(2) - 1)*(3*xi(2) - 2))/2;
      case 3
        vBase = (xi(3)*(3*xi(3) - 1)*(3*xi(3) - 2))/2;
      case 4
        vBase = (xi(4)*(3*xi(4) - 1)*(3*xi(4) - 2))/2;
      case 5
        vBase = (9*xi(1)*xi(2)*(3*xi(1) - 1))/2;
      case 6
        vBase = (9*xi(2)*xi(3)*(3*xi(2) - 1))/2;
      case 7
        vBase = (9*xi(1)*xi(3)*(3*xi(3) - 1))/2;
      case 8
        vBase = (9*xi(1)*xi(4)*(3*xi(1) - 1))/2;
      case 9
        vBase = (9*xi(2)*xi(4)*(3*xi(2) - 1))/2;
      case 10
        vBase = (9*xi(3)*xi(4)*(3*xi(3) - 1))/2;
      case 11
        vBase = (9*xi(1)*xi(2)*(3*xi(2) - 1))/2;
      case 12
        vBase = (9*xi(2)*xi(3)*(3*xi(3) - 1))/2;
      case 13
        vBase = (9*xi(1)*xi(3)*(3*xi(1) - 1))/2;
      case 14
        vBase = (9*xi(1)*xi(4)*(3*xi(4) - 1))/2;
      case 15
        vBase = (9*xi(2)*xi(4)*(3*xi(4) - 1))/2;
      case 16
        vBase = (9*xi(3)*xi(4)*(3*xi(4) - 1))/2;
      case 17
        vBase = 27*xi(1)*xi(2)*xi(3);
      case 18
        vBase = 27*xi(1)*xi(2)*xi(4);
      case 19
        vBase = 27*xi(2)*xi(3)*xi(4);
      case 20
        vBase = 27*xi(1)*xi(3)*xi(4);
    end

  case {2,3,4}

    switch i_dof
      case 1
        dNdxi1 = (27*xi(1)^2)/2 - 9*xi(1) + 1;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 2
        dNdxi1 = 0;
        dNdxi2 = (27*xi(2)^2)/2 - 9*xi(2) + 1;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 3
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (27*xi(3)^2)/2 - 9*xi(3) + 1;
        dNdxi4 = 0;
      case 4
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (27*xi(4)^2)/2 - 9*xi(4) + 1;
      case 5
        dNdxi1 = (9*xi(2)*(6*xi(1) - 1))/2;
        dNdxi2 = (9*xi(1)*(3*xi(1) - 1))/2;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 6
        dNdxi1 = 0;
        dNdxi2 = (9*xi(3)*(6*xi(2) - 1))/2;
        dNdxi3 = (9*xi(2)*(3*xi(2) - 1))/2;
        dNdxi4 = 0;
      case 7
        dNdxi1 = (9*xi(3)*(3*xi(3) - 1))/2;
        dNdxi2 = 0;
        dNdxi3 = (9*xi(1)*(6*xi(3) - 1))/2;
        dNdxi4 = 0;
      case 8
        dNdxi1 = (9*xi(4)*(6*xi(1) - 1))/2;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (9*xi(1)*(3*xi(1) - 1))/2;
      case 9
        dNdxi1 = 0;
        dNdxi2 = (9*xi(4)*(6*xi(2) - 1))/2;
        dNdxi3 = 0;
        dNdxi4 = (9*xi(2)*(3*xi(2) - 1))/2;
      case 10
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (9*xi(4)*(6*xi(3) - 1))/2;
        dNdxi4 = (9*xi(3)*(3*xi(3) - 1))/2;
      case 11
        dNdxi1 = (9*xi(2)*(3*xi(2) - 1))/2;
        dNdxi2 = (9*xi(1)*(6*xi(2) - 1))/2;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 12
        dNdxi1 = 0;
        dNdxi2 = (9*xi(3)*(3*xi(3) - 1))/2;
        dNdxi3 = (9*xi(2)*(6*xi(3) - 1))/2;
        dNdxi4 = 0;
      case 13
        dNdxi1 = (9*xi(3)*(6*xi(1) - 1))/2;
        dNdxi2 = 0;
        dNdxi3 = (9*xi(1)*(3*xi(1) - 1))/2;
        dNdxi4 = 0;
      case 14
        dNdxi1 = (9*xi(4)*(3*xi(4) - 1))/2;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (9*xi(1)*(6*xi(4) - 1))/2;
      case 15
        dNdxi1 = 0;
        dNdxi2 = (9*xi(4)*(3*xi(4) - 1))/2;
        dNdxi3 = 0;
        dNdxi4 = (9*xi(2)*(6*xi(4) - 1))/2;
      case 16
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (9*xi(4)*(3*xi(4) - 1))/2;
        dNdxi4 = (9*xi(3)*(6*xi(4) - 1))/2;
      case 17
        dNdxi1 = 27*xi(2)*xi(3);
        dNdxi2 = 27*xi(1)*xi(3);
        dNdxi3 = 27*xi(1)*xi(2);
        dNdxi4 = 0;
      case 18
        dNdxi1 = 27*xi(2)*xi(4);
        dNdxi2 = 27*xi(1)*xi(4);
        dNdxi3 = 0;
        dNdxi4 = 27*xi(1)*xi(2);
      case 19
        dNdxi1 = 0;
        dNdxi2 = 27*xi(3)*xi(4);
        dNdxi3 = 27*xi(2)*xi(4);
        dNdxi4 = 27*xi(2)*xi(3);
      case 20
        dNdxi1 = 27*xi(3)*xi(4);
        dNdxi2 = 0;
        dNdxi3 = 27*xi(1)*xi(4);
        dNdxi4 = 27*xi(1)*xi(3);
    end

    if( i_eval==2 )

      vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:,2)*dNdxi2 + aInvJac(:,3)*dNdxi3 + aInvJac(:,4)*dNdxi4;

    elseif( i_eval==3 )

      vBase = aInvJac(:,5)*dNdxi1 + aInvJac(:,6)*dNdxi2 + aInvJac(:,7)*dNdxi3 + aInvJac(:,8)*dNdxi4;

    else

      vBase = aInvJac(:,9)*dNdxi1 + aInvJac(:,10)*dNdxi2 + aInvJac(:,11)*dNdxi3 + aInvJac(:,12)*dNdxi4;

    end

  case {22,23,24,32,33,34,42,43,44}   % Evaluation of second derivatives.

    switch i_dof

      case 1
        d2Ndxi1dxi1 = 27*xi(1) - 9;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 2
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 27*xi(2) - 9;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 3
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 27*xi(3) - 9;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 4
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 27*xi(4) - 9;

      case 5
        d2Ndxi1dxi1 = 27*xi(2);
        d2Ndxi2dxi1 = 27*xi(1) - 9/2;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 27*xi(1) - 9/2;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 6
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 27*xi(3);
        d2Ndxi3dxi2 = 27*xi(2) - 9/2;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 27*xi(2) - 9/2;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 7
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 27*xi(3) - 9/2;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 27*xi(3) - 9/2;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 27*xi(1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 8
        d2Ndxi1dxi1 = 27*xi(4);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 27*xi(1) - 9/2;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 27*xi(1) - 9/2;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 9
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 27*xi(4);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 27*xi(2) - 9/2;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 27*xi(2) - 9/2;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 10
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 27*xi(4);
        d2Ndxi4dxi3 = 27*xi(3) - 9/2;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 27*xi(3) - 9/2;
        d2Ndxi4dxi4 = 0;

      case 11
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 27*xi(2) - 9/2;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 27*xi(2) - 9/2;
        d2Ndxi2dxi2 = 27*xi(1);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 12
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 27*xi(3) - 9/2;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 27*xi(3) - 9/2;
        d2Ndxi3dxi3 = 27*xi(2);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 13
        d2Ndxi1dxi1 = 27*xi(3);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 27*xi(1) - 9/2;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 27*xi(1) - 9/2;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 14
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 27*xi(4) - 9/2;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 27*xi(4) - 9/2;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 27*xi(1);

      case 15
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 27*xi(4) - 9/2;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 27*xi(4) - 9/2;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 27*xi(2);

      case 16
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 27*xi(4) - 9/2;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 27*xi(4) - 9/2;
        d2Ndxi4dxi4 = 27*xi(3);

      case 17
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 27*xi(3);
        d2Ndxi3dxi1 = 27*xi(2);
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 27*xi(3);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 27*xi(1);
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 27*xi(2);
        d2Ndxi2dxi3 = 27*xi(1);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 18
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 27*xi(4);
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 27*xi(2);
        d2Ndxi1dxi2 = 27*xi(4);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 27*xi(1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 27*xi(2);
        d2Ndxi2dxi4 = 27*xi(1);
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 19
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 27*xi(4);
        d2Ndxi4dxi2 = 27*xi(3);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 27*xi(4);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 27*xi(2);
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 27*xi(3);
        d2Ndxi3dxi4 = 27*xi(2);
        d2Ndxi4dxi4 = 0;

      case 20
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 27*xi(4);
        d2Ndxi4dxi1 = 27*xi(3);
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 27*xi(4);
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 27*xi(1);
        d2Ndxi1dxi4 = 27*xi(3);
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 27*xi(1);
        d2Ndxi4dxi4 = 0;
    end

    switch( i_eval )
      case 22
        vBase = aInvJac(:,1).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:, 2)*d2Ndxi2dxi1 + aInvJac(:, 3)*d2Ndxi3dxi1 + aInvJac(:, 4)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,2).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:, 2)*d2Ndxi2dxi2 + aInvJac(:, 3)*d2Ndxi3dxi2 + aInvJac(:, 4)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,3).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:, 2)*d2Ndxi2dxi3 + aInvJac(:, 3)*d2Ndxi3dxi3 + aInvJac(:, 4)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,4).*( aInvJac(:,1)*d2Ndxi1dxi4 + aInvJac(:, 2)*d2Ndxi2dxi4 + aInvJac(:, 3)*d2Ndxi3dxi4 + aInvJac(:, 4)*d2Ndxi4dxi4 );
      case 33
        vBase = aInvJac(:,5).*( aInvJac(:,5)*d2Ndxi1dxi1 + aInvJac(:, 6)*d2Ndxi2dxi1 + aInvJac(:, 7)*d2Ndxi3dxi1 + aInvJac(:, 8)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,6).*( aInvJac(:,5)*d2Ndxi1dxi2 + aInvJac(:, 6)*d2Ndxi2dxi2 + aInvJac(:, 7)*d2Ndxi3dxi2 + aInvJac(:, 8)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,7).*( aInvJac(:,5)*d2Ndxi1dxi3 + aInvJac(:, 6)*d2Ndxi2dxi3 + aInvJac(:, 7)*d2Ndxi3dxi3 + aInvJac(:, 8)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,8).*( aInvJac(:,5)*d2Ndxi1dxi4 + aInvJac(:, 6)*d2Ndxi2dxi4 + aInvJac(:, 7)*d2Ndxi3dxi4 + aInvJac(:, 8)*d2Ndxi4dxi4 );

      case 44
        vBase = aInvJac(:, 9).*( aInvJac(:,9)*d2Ndxi1dxi1 + aInvJac(:,10)*d2Ndxi2dxi1 + aInvJac(:,11)*d2Ndxi3dxi1 + aInvJac(:,12)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,10).*( aInvJac(:,9)*d2Ndxi1dxi2 + aInvJac(:,10)*d2Ndxi2dxi2 + aInvJac(:,11)*d2Ndxi3dxi2 + aInvJac(:,12)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,11).*( aInvJac(:,9)*d2Ndxi1dxi3 + aInvJac(:,10)*d2Ndxi2dxi3 + aInvJac(:,11)*d2Ndxi3dxi3 + aInvJac(:,12)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,12).*( aInvJac(:,9)*d2Ndxi1dxi4 + aInvJac(:,10)*d2Ndxi2dxi4 + aInvJac(:,11)*d2Ndxi3dxi4 + aInvJac(:,12)*d2Ndxi4dxi4 );

      case {23,32}
        vBase = aInvJac(:,5).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:, 2)*d2Ndxi2dxi1 + aInvJac(:, 3)*d2Ndxi3dxi1 + aInvJac(:, 4)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,6).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:, 2)*d2Ndxi2dxi2 + aInvJac(:, 3)*d2Ndxi3dxi2 + aInvJac(:, 4)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,7).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:, 2)*d2Ndxi2dxi3 + aInvJac(:, 3)*d2Ndxi3dxi3 + aInvJac(:, 4)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,8).*( aInvJac(:,1)*d2Ndxi1dxi4 + aInvJac(:, 2)*d2Ndxi2dxi4 + aInvJac(:, 3)*d2Ndxi3dxi4 + aInvJac(:, 4)*d2Ndxi4dxi4 );

      case {24,42}
        vBase = aInvJac(:, 9).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:, 2)*d2Ndxi2dxi1 + aInvJac(:, 3)*d2Ndxi3dxi1 + aInvJac(:, 4)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,10).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:, 2)*d2Ndxi2dxi2 + aInvJac(:, 3)*d2Ndxi3dxi2 + aInvJac(:, 4)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,11).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:, 2)*d2Ndxi2dxi3 + aInvJac(:, 3)*d2Ndxi3dxi3 + aInvJac(:, 4)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,12).*( aInvJac(:,1)*d2Ndxi1dxi4 + aInvJac(:, 2)*d2Ndxi2dxi4 + aInvJac(:, 3)*d2Ndxi3dxi4 + aInvJac(:, 4)*d2Ndxi4dxi4 );

      case {34,43}
        vBase = aInvJac(:, 9).*( aInvJac(:,5)*d2Ndxi1dxi1 + aInvJac(:, 6)*d2Ndxi2dxi1 + aInvJac(:, 7)*d2Ndxi3dxi1 + aInvJac(:, 8)*d2Ndxi4dxi1 ) + ...
                aInvJac(:,10).*( aInvJac(:,5)*d2Ndxi1dxi2 + aInvJac(:, 6)*d2Ndxi2dxi2 + aInvJac(:, 7)*d2Ndxi3dxi2 + aInvJac(:, 8)*d2Ndxi4dxi2 ) + ...
                aInvJac(:,11).*( aInvJac(:,5)*d2Ndxi1dxi3 + aInvJac(:, 6)*d2Ndxi2dxi3 + aInvJac(:, 7)*d2Ndxi3dxi3 + aInvJac(:, 8)*d2Ndxi4dxi3 ) + ...
                aInvJac(:,12).*( aInvJac(:,5)*d2Ndxi1dxi4 + aInvJac(:, 6)*d2Ndxi2dxi4 + aInvJac(:, 7)*d2Ndxi3dxi4 + aInvJac(:, 8)*d2Ndxi4dxi4 );
    end

  otherwise

    vBase = 0;

end
