function [ vBase, nLDof, xLDof, sfun ] = sf_tet_P4( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_TET_P4 Fourth order Lagrange shape functions for tetrahedra (P4).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_TET_P4( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fourth order P4 Lagrange shape functions on 3D tetrahedral elements
%   with values defined in the nodes, edges, faces, and cell center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 3              Number of space dimensions
%       n_vert      scalar: 4              Number of vertices per cell
%       i_dof       scalar: 1-35           Local basis function to evaluate
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
%   See also SF_TET_P1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [4 18 12 1];
xLDof = [1 0 0 0 3/4   0 1/4 3/4   0   0 1/2   0 1/2 1/2   0   0 1/4   0 3/4 1/4   0   0 1/2 1/2   0 1/4 1/4 1/4   0 1/2 1/4 1/4   0 1/4 1/4;
         0 1 0 0 1/4 3/4   0   0 3/4   0 1/2 1/2   0   0 1/2   0 3/4 1/4   0   0 1/4   0 1/4 1/4 1/2   0 1/2 1/2 1/4   0 1/4 1/4 1/4   0 1/4;
         0 0 1 0   0 1/4 3/4   0   0 3/4   0 1/2 1/2   0   0 1/2   0 3/4 1/4   0   0 1/4 1/4   0 1/4 1/2 1/4   0 1/2 1/4 1/2   0 1/4 1/4 1/4;
         0 0 0 1   0   0   0 1/4 1/4 1/4   0   0   0 1/2 1/2 1/2   0   0   0 3/4 3/4 3/4   0 1/4 1/4 1/4   0 1/4 1/4 1/4   0 1/2 1/2 1/2 1/4];
sfun  = 'sf_tet_P4';


switch i_eval

  case 1

    switch i_dof
      case 1
        vBase = (xi(1)*(2*xi(1) - 1)*(4*xi(1) - 1)*(4*xi(1) - 3))/3;
      case 2
        vBase = (xi(2)*(2*xi(2) - 1)*(4*xi(2) - 1)*(4*xi(2) - 3))/3;
      case 3
        vBase = (xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1)*(4*xi(3) - 3))/3;
      case 4
        vBase = (xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1)*(4*xi(4) - 3))/3;
      case 5
        vBase = (16*xi(1)*xi(2)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
      case 6
        vBase = (16*xi(2)*xi(3)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
      case 7
        vBase = (16*xi(1)*xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
      case 8
        vBase = (16*xi(1)*xi(4)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
      case 9
        vBase = (16*xi(2)*xi(4)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
      case 10
        vBase = (16*xi(3)*xi(4)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
      case 11
        vBase = 4*xi(1)*xi(2)*(4*xi(1) - 1)*(4*xi(2) - 1);
      case 12
        vBase = 4*xi(2)*xi(3)*(4*xi(2) - 1)*(4*xi(3) - 1);
      case 13
        vBase = 4*xi(1)*xi(3)*(4*xi(1) - 1)*(4*xi(3) - 1);
      case 14
        vBase = 4*xi(1)*xi(4)*(4*xi(1) - 1)*(4*xi(4) - 1);
      case 15
        vBase = 4*xi(2)*xi(4)*(4*xi(2) - 1)*(4*xi(4) - 1);
      case 16
        vBase = 4*xi(3)*xi(4)*(4*xi(3) - 1)*(4*xi(4) - 1);
      case 17
        vBase = (16*xi(1)*xi(2)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
      case 18
        vBase = (16*xi(2)*xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
      case 19
        vBase = (16*xi(1)*xi(3)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
      case 20
        vBase = (16*xi(1)*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
      case 21
        vBase = (16*xi(2)*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
      case 22
        vBase = (16*xi(3)*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
      case 23
        vBase = 32*xi(1)*xi(2)*xi(3)*(4*xi(1) - 1);
      case 24
        vBase = 32*xi(1)*xi(2)*xi(4)*(4*xi(1) - 1);
      case 25
        vBase = 32*xi(2)*xi(3)*xi(4)*(4*xi(2) - 1);
      case 26
        vBase = 32*xi(1)*xi(3)*xi(4)*(4*xi(3) - 1);
      case 27
        vBase = 32*xi(1)*xi(2)*xi(3)*(4*xi(2) - 1);
      case 28
        vBase = 32*xi(1)*xi(2)*xi(4)*(4*xi(2) - 1);
      case 29
        vBase = 32*xi(2)*xi(3)*xi(4)*(4*xi(3) - 1);
      case 30
        vBase = 32*xi(1)*xi(3)*xi(4)*(4*xi(1) - 1);
      case 31
        vBase = 32*xi(1)*xi(2)*xi(3)*(4*xi(3) - 1);
      case 32
        vBase = 32*xi(1)*xi(2)*xi(4)*(4*xi(4) - 1);
      case 33
        vBase = 32*xi(2)*xi(3)*xi(4)*(4*xi(4) - 1);
      case 34
        vBase = 32*xi(1)*xi(3)*xi(4)*(4*xi(4) - 1);
      case 35
        vBase = 256*xi(1)*xi(2)*xi(3)*xi(4);
    end

  case {2,3,4}

    switch i_dof

      case 1
        dNdxi1 = ((8*xi(1) - 3)*(16*xi(1)^2 - 12*xi(1) + 1))/3;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 2
        dNdxi1 = 0;
        dNdxi2 = ((8*xi(2) - 3)*(16*xi(2)^2 - 12*xi(2) + 1))/3;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 3
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = ((8*xi(3) - 3)*(16*xi(3)^2 - 12*xi(3) + 1))/3;
        dNdxi4 = 0;
      case 4
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = ((8*xi(4) - 3)*(16*xi(4)^2 - 12*xi(4) + 1))/3;
      case 5
        dNdxi1 = (16*xi(2)*(24*xi(1)^2 - 12*xi(1) + 1))/3;
        dNdxi2 = (16*xi(1)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 6
        dNdxi1 = 0;
        dNdxi2 = (16*xi(3)*(24*xi(2)^2 - 12*xi(2) + 1))/3;
        dNdxi3 = (16*xi(2)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
        dNdxi4 = 0;
      case 7
        dNdxi1 = (16*xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
        dNdxi2 = 0;
        dNdxi3 = (16*xi(1)*(24*xi(3)^2 - 12*xi(3) + 1))/3;
        dNdxi4 = 0;
      case 8
        dNdxi1 = (16*xi(4)*(24*xi(1)^2 - 12*xi(1) + 1))/3;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (16*xi(1)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
      case 9
        dNdxi1 = 0;
        dNdxi2 = (16*xi(4)*(24*xi(2)^2 - 12*xi(2) + 1))/3;
        dNdxi3 = 0;
        dNdxi4 = (16*xi(2)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
      case 10
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (16*xi(4)*(24*xi(3)^2 - 12*xi(3) + 1))/3;
        dNdxi4 = (16*xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
      case 11
        dNdxi1 = 4*xi(2)*(8*xi(1) - 1)*(4*xi(2) - 1);
        dNdxi2 = 4*xi(1)*(4*xi(1) - 1)*(8*xi(2) - 1);
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 12
        dNdxi1 = 0;
        dNdxi2 = 4*xi(3)*(8*xi(2) - 1)*(4*xi(3) - 1);
        dNdxi3 = 4*xi(2)*(4*xi(2) - 1)*(8*xi(3) - 1);
        dNdxi4 = 0;
      case 13
        dNdxi1 = 4*xi(3)*(8*xi(1) - 1)*(4*xi(3) - 1);
        dNdxi2 = 0;
        dNdxi3 = 4*xi(1)*(4*xi(1) - 1)*(8*xi(3) - 1);
        dNdxi4 = 0;
      case 14
        dNdxi1 = 4*xi(4)*(8*xi(1) - 1)*(4*xi(4) - 1);
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = 4*xi(1)*(4*xi(1) - 1)*(8*xi(4) - 1);
      case 15
        dNdxi1 = 0;
        dNdxi2 = 4*xi(4)*(8*xi(2) - 1)*(4*xi(4) - 1);
        dNdxi3 = 0;
        dNdxi4 = 4*xi(2)*(4*xi(2) - 1)*(8*xi(4) - 1);
      case 16
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = 4*xi(4)*(8*xi(3) - 1)*(4*xi(4) - 1);
        dNdxi4 = 4*xi(3)*(4*xi(3) - 1)*(8*xi(4) - 1);
      case 17
        dNdxi1 = (16*xi(2)*(2*xi(2) - 1)*(4*xi(2) - 1))/3;
        dNdxi2 = (16*xi(1)*(24*xi(2)^2 - 12*xi(2) + 1))/3;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 18
        dNdxi1 = 0;
        dNdxi2 = (16*xi(3)*(2*xi(3) - 1)*(4*xi(3) - 1))/3;
        dNdxi3 = (16*xi(2)*(24*xi(3)^2 - 12*xi(3) + 1))/3;
        dNdxi4 = 0;
      case 19
        dNdxi1 = (16*xi(3)*(24*xi(1)^2 - 12*xi(1) + 1))/3;
        dNdxi2 = 0;
        dNdxi3 = (16*xi(1)*(2*xi(1) - 1)*(4*xi(1) - 1))/3;
        dNdxi4 = 0;
      case 20
        dNdxi1 = (16*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (16*xi(1)*(24*xi(4)^2 - 12*xi(4) + 1))/3;
      case 21
        dNdxi1 = 0;
        dNdxi2 = (16*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
        dNdxi3 = 0;
        dNdxi4 = (16*xi(2)*(24*xi(4)^2 - 12*xi(4) + 1))/3;
      case 22
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (16*xi(4)*(2*xi(4) - 1)*(4*xi(4) - 1))/3;
        dNdxi4 = (16*xi(3)*(24*xi(4)^2 - 12*xi(4) + 1))/3;
      case 23
        dNdxi1 = 32*xi(2)*xi(3)*(8*xi(1) - 1);
        dNdxi2 = 32*xi(1)*xi(3)*(4*xi(1) - 1);
        dNdxi3 = 32*xi(1)*xi(2)*(4*xi(1) - 1);
        dNdxi4 = 0;
      case 24
        dNdxi1 = 32*xi(2)*xi(4)*(8*xi(1) - 1);
        dNdxi2 = 32*xi(1)*xi(4)*(4*xi(1) - 1);
        dNdxi3 = 0;
        dNdxi4 = 32*xi(1)*xi(2)*(4*xi(1) - 1);
      case 25
        dNdxi1 = 0;
        dNdxi2 = 32*xi(3)*xi(4)*(8*xi(2) - 1);
        dNdxi3 = 32*xi(2)*xi(4)*(4*xi(2) - 1);
        dNdxi4 = 32*xi(2)*xi(3)*(4*xi(2) - 1);
      case 26
        dNdxi1 = 32*xi(3)*xi(4)*(4*xi(3) - 1);
        dNdxi2 = 0;
        dNdxi3 = 32*xi(1)*xi(4)*(8*xi(3) - 1);
        dNdxi4 = 32*xi(1)*xi(3)*(4*xi(3) - 1);
      case 27
        dNdxi1 = 32*xi(2)*xi(3)*(4*xi(2) - 1);
        dNdxi2 = 32*xi(1)*xi(3)*(8*xi(2) - 1);
        dNdxi3 = 32*xi(1)*xi(2)*(4*xi(2) - 1);
        dNdxi4 = 0;
      case 28
        dNdxi1 = 32*xi(2)*xi(4)*(4*xi(2) - 1);
        dNdxi2 = 32*xi(1)*xi(4)*(8*xi(2) - 1);
        dNdxi3 = 0;
        dNdxi4 = 32*xi(1)*xi(2)*(4*xi(2) - 1);
      case 29
        dNdxi1 = 0;
        dNdxi2 = 32*xi(3)*xi(4)*(4*xi(3) - 1);
        dNdxi3 = 32*xi(2)*xi(4)*(8*xi(3) - 1);
        dNdxi4 = 32*xi(2)*xi(3)*(4*xi(3) - 1);
      case 30
        dNdxi1 = 32*xi(3)*xi(4)*(8*xi(1) - 1);
        dNdxi2 = 0;
        dNdxi3 = 32*xi(1)*xi(4)*(4*xi(1) - 1);
        dNdxi4 = 32*xi(1)*xi(3)*(4*xi(1) - 1);
      case 31
        dNdxi1 = 32*xi(2)*xi(3)*(4*xi(3) - 1);
        dNdxi2 = 32*xi(1)*xi(3)*(4*xi(3) - 1);
        dNdxi3 = 32*xi(1)*xi(2)*(8*xi(3) - 1);
        dNdxi4 = 0;
      case 32
        dNdxi1 = 32*xi(2)*xi(4)*(4*xi(4) - 1);
        dNdxi2 = 32*xi(1)*xi(4)*(4*xi(4) - 1);
        dNdxi3 = 0;
        dNdxi4 = 32*xi(1)*xi(2)*(8*xi(4) - 1);
      case 33
        dNdxi1 = 0;
        dNdxi2 = 32*xi(3)*xi(4)*(4*xi(4) - 1);
        dNdxi3 = 32*xi(2)*xi(4)*(4*xi(4) - 1);
        dNdxi4 = 32*xi(2)*xi(3)*(8*xi(4) - 1);
      case 34
        dNdxi1 = 32*xi(3)*xi(4)*(4*xi(4) - 1);
        dNdxi2 = 0;
        dNdxi3 = 32*xi(1)*xi(4)*(4*xi(4) - 1);
        dNdxi4 = 32*xi(1)*xi(3)*(8*xi(4) - 1);
      case 35
        dNdxi1 = 256*xi(2)*xi(3)*xi(4);
        dNdxi2 = 256*xi(1)*xi(3)*xi(4);
        dNdxi3 = 256*xi(1)*xi(2)*xi(4);
        dNdxi4 = 256*xi(1)*xi(2)*xi(3);
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
        d2Ndxi1dxi1 = 128*xi(1)^2 - 96*xi(1) + 44/3;
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
        d2Ndxi2dxi2 = 128*xi(2)^2 - 96*xi(2) + 44/3;
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
        d2Ndxi3dxi3 = 128*xi(3)^2 - 96*xi(3) + 44/3;
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
        d2Ndxi4dxi4 = 128*xi(4)^2 - 96*xi(4) + 44/3;

      case 5
        d2Ndxi1dxi1 = 64*xi(2)*(4*xi(1) - 1);
        d2Ndxi2dxi1 = 128*xi(1)^2 - 64*xi(1) + 16/3;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 128*xi(1)^2 - 64*xi(1) + 16/3;
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
        d2Ndxi2dxi2 = 64*xi(3)*(4*xi(2) - 1);
        d2Ndxi3dxi2 = 128*xi(2)^2 - 64*xi(2) + 16/3;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 128*xi(2)^2 - 64*xi(2) + 16/3;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 7
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 64*xi(1)*(4*xi(3) - 1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 8
        d2Ndxi1dxi1 = 64*xi(4)*(4*xi(1) - 1);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 128*xi(1)^2 - 64*xi(1) + 16/3;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 128*xi(1)^2 - 64*xi(1) + 16/3;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 9
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 64*xi(4)*(4*xi(2) - 1);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 128*xi(2)^2 - 64*xi(2) + 16/3;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 128*xi(2)^2 - 64*xi(2) + 16/3;
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
        d2Ndxi3dxi3 = 64*xi(4)*(4*xi(3) - 1);
        d2Ndxi4dxi3 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi4dxi4 = 0;

      case 11
        d2Ndxi1dxi1 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi2dxi1 = 4*(8*xi(1) - 1)*(8*xi(2) - 1);
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 4*(8*xi(1) - 1)*(8*xi(2) - 1);
        d2Ndxi2dxi2 = 32*xi(1)*(4*xi(1) - 1);
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
        d2Ndxi2dxi2 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi3dxi2 = 4*(8*xi(2) - 1)*(8*xi(3) - 1);
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 4*(8*xi(2) - 1)*(8*xi(3) - 1);
        d2Ndxi3dxi3 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 13
        d2Ndxi1dxi1 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 4*(8*xi(1) - 1)*(8*xi(3) - 1);
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 4*(8*xi(1) - 1)*(8*xi(3) - 1);
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 14
        d2Ndxi1dxi1 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 4*(8*xi(1) - 1)*(8*xi(4) - 1);
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 4*(8*xi(1) - 1)*(8*xi(4) - 1);
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 32*xi(1)*(4*xi(1) - 1);

      case 15
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 4*(8*xi(2) - 1)*(8*xi(4) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 4*(8*xi(2) - 1)*(8*xi(4) - 1);
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 32*xi(2)*(4*xi(2) - 1);

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
        d2Ndxi3dxi3 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi4dxi3 = 4*(8*xi(3) - 1)*(8*xi(4) - 1);
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 4*(8*xi(3) - 1)*(8*xi(4) - 1);
        d2Ndxi4dxi4 = 32*xi(3)*(4*xi(3) - 1);

      case 17
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 128*xi(2)^2 - 64*xi(2) + 16/3;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 128*xi(2)^2 - 64*xi(2) + 16/3;
        d2Ndxi2dxi2 = 64*xi(1)*(4*xi(2) - 1);
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

      case 18
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 128*xi(3)^2 - 64*xi(3) + 16/3;
        d2Ndxi3dxi3 = 64*xi(2)*(4*xi(3) - 1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 19
        d2Ndxi1dxi1 = 64*xi(3)*(4*xi(1) - 1);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 128*xi(1)^2 - 64*xi(1) + 16/3;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 128*xi(1)^2 - 64*xi(1) + 16/3;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 20
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 64*xi(1)*(4*xi(4) - 1);

      case 21
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 64*xi(2)*(4*xi(4) - 1);

      case 22
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
        d2Ndxi4dxi3 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 128*xi(4)^2 - 64*xi(4) + 16/3;
        d2Ndxi4dxi4 = 64*xi(3)*(4*xi(4) - 1);

      case 23
        d2Ndxi1dxi1 = 256*xi(2)*xi(3);
        d2Ndxi2dxi1 = 32*xi(3)*(8*xi(1) - 1);
        d2Ndxi3dxi1 = 32*xi(2)*(8*xi(1) - 1);
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 32*xi(3)*(8*xi(1) - 1);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(2)*(8*xi(1) - 1);
        d2Ndxi2dxi3 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 24
        d2Ndxi1dxi1 = 256*xi(2)*xi(4);
        d2Ndxi2dxi1 = 32*xi(4)*(8*xi(1) - 1);
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 32*xi(2)*(8*xi(1) - 1);
        d2Ndxi1dxi2 = 32*xi(4)*(8*xi(1) - 1);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 32*xi(2)*(8*xi(1) - 1);
        d2Ndxi2dxi4 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 25
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 256*xi(3)*xi(4);
        d2Ndxi3dxi2 = 32*xi(4)*(8*xi(2) - 1);
        d2Ndxi4dxi2 = 32*xi(3)*(8*xi(2) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 32*xi(4)*(8*xi(2) - 1);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 32*xi(3)*(8*xi(2) - 1);
        d2Ndxi3dxi4 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi4dxi4 = 0;

      case 26
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 32*xi(4)*(8*xi(3) - 1);
        d2Ndxi4dxi1 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(4)*(8*xi(3) - 1);
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 256*xi(1)*xi(4);
        d2Ndxi4dxi3 = 32*xi(1)*(8*xi(3) - 1);
        d2Ndxi1dxi4 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 32*xi(1)*(8*xi(3) - 1);
        d2Ndxi4dxi4 = 0;

      case 27
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 32*xi(3)*(8*xi(2) - 1);
        d2Ndxi3dxi1 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 32*xi(3)*(8*xi(2) - 1);
        d2Ndxi2dxi2 = 256*xi(1)*xi(3);
        d2Ndxi3dxi2 = 32*xi(1)*(8*xi(2) - 1);
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi2dxi3 = 32*xi(1)*(8*xi(2) - 1);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 28
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 32*xi(4)*(8*xi(2) - 1);
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi1dxi2 = 32*xi(4)*(8*xi(2) - 1);
        d2Ndxi2dxi2 = 256*xi(1)*xi(4);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 32*xi(1)*(8*xi(2) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 32*xi(2)*(4*xi(2) - 1);
        d2Ndxi2dxi4 = 32*xi(1)*(8*xi(2) - 1);
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 29
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 32*xi(4)*(8*xi(3) - 1);
        d2Ndxi4dxi2 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 32*xi(4)*(8*xi(3) - 1);
        d2Ndxi3dxi3 = 256*xi(2)*xi(4);
        d2Ndxi4dxi3 = 32*xi(2)*(8*xi(3) - 1);
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi3dxi4 = 32*xi(2)*(8*xi(3) - 1);
        d2Ndxi4dxi4 = 0;

      case 30
        d2Ndxi1dxi1 = 256*xi(3)*xi(4);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 32*xi(4)*(8*xi(1) - 1);
        d2Ndxi4dxi1 = 32*xi(3)*(8*xi(1) - 1);
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(4)*(8*xi(1) - 1);
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi1dxi4 = 32*xi(3)*(8*xi(1) - 1);
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 32*xi(1)*(4*xi(1) - 1);
        d2Ndxi4dxi4 = 0;

      case 31
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi3dxi1 = 32*xi(2)*(8*xi(3) - 1);
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 32*xi(3)*(4*xi(3) - 1);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 32*xi(1)*(8*xi(3) - 1);
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(2)*(8*xi(3) - 1);
        d2Ndxi2dxi3 = 32*xi(1)*(8*xi(3) - 1);
        d2Ndxi3dxi3 = 256*xi(1)*xi(2);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 32
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 32*xi(2)*(8*xi(4) - 1);
        d2Ndxi1dxi2 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 32*xi(1)*(8*xi(4) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 32*xi(2)*(8*xi(4) - 1);
        d2Ndxi2dxi4 = 32*xi(1)*(8*xi(4) - 1);
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 256*xi(1)*xi(2);

      case 33
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi4dxi2 = 32*xi(3)*(8*xi(4) - 1);
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 32*xi(2)*(8*xi(4) - 1);
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 32*xi(3)*(8*xi(4) - 1);
        d2Ndxi3dxi4 = 32*xi(2)*(8*xi(4) - 1);
        d2Ndxi4dxi4 = 256*xi(2)*xi(3);

      case 34
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi4dxi1 = 32*xi(3)*(8*xi(4) - 1);
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 32*xi(4)*(4*xi(4) - 1);
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 32*xi(1)*(8*xi(4) - 1);
        d2Ndxi1dxi4 = 32*xi(3)*(8*xi(4) - 1);
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 32*xi(1)*(8*xi(4) - 1);
        d2Ndxi4dxi4 = 256*xi(1)*xi(3);

      case 35
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 256*xi(3)*xi(4);
        d2Ndxi3dxi1 = 256*xi(2)*xi(4);
        d2Ndxi4dxi1 = 256*xi(2)*xi(3);
        d2Ndxi1dxi2 = 256*xi(3)*xi(4);
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 256*xi(1)*xi(4);
        d2Ndxi4dxi2 = 256*xi(1)*xi(3);
        d2Ndxi1dxi3 = 256*xi(2)*xi(4);
        d2Ndxi2dxi3 = 256*xi(1)*xi(4);
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 256*xi(1)*xi(2);
        d2Ndxi1dxi4 = 256*xi(2)*xi(3);
        d2Ndxi2dxi4 = 256*xi(1)*xi(3);
        d2Ndxi3dxi4 = 256*xi(1)*xi(2);
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
