function [ vBase, nLDof, xLDof, sfun ] = sf_tri_P5( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_TRI_P5 Fifth order Lagrange shape functions for triangles (P5).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_TRI_P5( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fifth order P5 Lagrange shape functions on 2D triangular elements
%   with values defined in the nodes, edges, and center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 2              Number of space dimensions
%       n_vert      scalar: 3              Number of vertices per cell
%       i_dof       scalar: 1-21           Local basis function to evaluate
%       xi          array [3,1]            Local coordinates of evaluation point
%       aInvJac     [n,6]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [3,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_TRI_P1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [3 12 0 6];
xLDof = [1 0 0 4/5   0 1/5 3/5   0 2/5 2/5   0 3/5 1/5   0 4/5 3/5 1/5 1/5 2/5 1/5 2/5;
         0 1 0 1/5 4/5   0 2/5 3/5   0 3/5 2/5   0 4/5 1/5   0 1/5 3/5 1/5 2/5 2/5 1/5;
         0 0 1   0 1/5 4/5   0 2/5 3/5   0 3/5 2/5   0 4/5 1/5 1/5 1/5 3/5 1/5 2/5 2/5];
sfun  = 'sf_tri_P5';


switch i_eval

  case 1

    switch i_dof
      case 1
        vBase = (xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3)*(5*xi(1) - 4))/24;
      case 2
        vBase = (xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3)*(5*xi(2) - 4))/24;
      case 3
        vBase = (xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3)*(5*xi(3) - 4))/24;
      case 4
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 5
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 6
        vBase = (25*xi(1)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 7
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(2) - 1))/12;
      case 8
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(3) - 1))/12;
      case 9
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
      case 10
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1)*(5*xi(2) - 2))/12;
      case 11
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
      case 12
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(3) - 1))/12;
      case 13
        vBase = (25*xi(1)*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 14
        vBase = (25*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 15
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 16
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 17
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 18
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 19
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
      case 20
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
      case 21
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
    end

  case {2,3}

    switch i_dof

      case 1
        dNdxi1 = (3125*xi(1)^4)/24 - (625*xi(1)^3)/3 + (875*xi(1)^2)/8 - (125*xi(1))/6 + 1;
        dNdxi2 = 0;
        dNdxi3 = 0;
      case 2
        dNdxi1 = 0;
        dNdxi2 = (3125*xi(2)^4)/24 - (625*xi(2)^3)/3 + (875*xi(2)^2)/8 - (125*xi(2))/6 + 1;
        dNdxi3 = 0;
      case 3
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (3125*xi(3)^4)/24 - (625*xi(3)^3)/3 + (875*xi(3)^2)/8 - (125*xi(3))/6 + 1;
      case 4
        dNdxi1 = (25*xi(2)*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
        dNdxi3 = 0;
      case 5
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 6
        dNdxi1 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
      case 7
        dNdxi1 = (25*xi(2)*(5*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(10*xi(2) - 1))/12;
        dNdxi3 = 0;
      case 8
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(5*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(10*xi(3) - 1))/12;
      case 9
        dNdxi1 = (25*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
      case 10
        dNdxi1 = (25*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1)*(5*xi(2) - 2))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        dNdxi3 = 0;
      case 11
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
      case 12
        dNdxi1 = (25*xi(3)*(5*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(10*xi(3) - 1))/12;
      case 13
        dNdxi1 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
        dNdxi2 = (25*xi(1)*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        dNdxi3 = 0;
      case 14
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
        dNdxi3 = (25*xi(2)*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
      case 15
        dNdxi1 = (25*xi(3)*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 16
        dNdxi1 = (125*xi(2)*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 17
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 18
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
      case 19
        dNdxi1 = (125*xi(2)*xi(3)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
      case 20
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
      case 21
        dNdxi1 = (125*xi(2)*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
    end

    if( i_eval==2 )

      vBase = aInvJac(:,1).*dNdxi1 + aInvJac(:,2).*dNdxi2 + aInvJac(:,3).*dNdxi3;

    else

      vBase = aInvJac(:,4).*dNdxi1 + aInvJac(:,5).*dNdxi2 + aInvJac(:,6).*dNdxi3;

    end

  case {22,23,32,33}   % Evaluation of second derivatives.

    switch i_dof

      case 1
        d2Ndxi1dxi1 = (125*(5*xi(1) - 2)*(10*xi(1)^2 - 8*xi(1) + 1))/12;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 2
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*(5*xi(2) - 2)*(10*xi(2)^2 - 8*xi(2) + 1))/12;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 3
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*(5*xi(3) - 2)*(10*xi(3)^2 - 8*xi(3) + 1))/12;

      case 4
        d2Ndxi1dxi1 = (125*xi(2)*(150*xi(1)^2 - 90*xi(1) + 11))/12;
        d2Ndxi2dxi1 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 5
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(3)*(150*xi(2)^2 - 90*xi(2) + 11))/12;
        d2Ndxi3dxi2 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi3dxi3 = 0;

      case 6
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(150*xi(3)^2 - 90*xi(3) + 11))/12;

      case 7
        d2Ndxi1dxi1 = (125*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi1 = (25*(10*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 8
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi2 = (25*(10*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;

      case 9
        d2Ndxi1dxi1 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(3) - 1))/2;

      case 10
        d2Ndxi1dxi1 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi2dxi1 = (25*(10*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(2) - 1))/2;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 11
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi3dxi2 = (25*(10*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(3) - 1))/2;

      case 12
        d2Ndxi1dxi1 = (125*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;

      case 13
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(150*xi(2)^2 - 90*xi(2) + 11))/12;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 14
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(150*xi(3)^2 - 90*xi(3) + 11))/12;

      case 15
        d2Ndxi1dxi1 = (125*xi(3)*(150*xi(1)^2 - 90*xi(1) + 11))/12;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;

      case 16
        d2Ndxi1dxi1 = 625*xi(2)*xi(3)*(5*xi(1) - 1);
        d2Ndxi2dxi1 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi1dxi2 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi1dxi3 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi3dxi3 = 0;

      case 17
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi1dxi2 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi2dxi2 = 625*xi(1)*xi(3)*(5*xi(2) - 1);
        d2Ndxi3dxi2 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi1dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi3 = 0;

      case 18
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi1dxi2 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi1dxi3 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi3dxi3 = 625*xi(1)*xi(2)*(5*xi(3) - 1);

      case 19
        d2Ndxi1dxi1 = (625*xi(2)*xi(3)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(3)*(5*xi(1) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi1dxi3 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi3 = 0;

      case 20
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi3 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(1)*xi(2)*(5*xi(2) - 1))/2;

      case 21
        d2Ndxi1dxi1 = (625*xi(2)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi3 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(1)*xi(2)*(5*xi(1) - 1))/2;
    end

    if( i_eval==22 )
      vBase = aInvJac(:,1).*( aInvJac(:,1).*d2Ndxi1dxi1 + aInvJac(:,2).*d2Ndxi2dxi1 + aInvJac(:,3).*d2Ndxi3dxi1 ) + ...
              aInvJac(:,2).*( aInvJac(:,1).*d2Ndxi1dxi2 + aInvJac(:,2).*d2Ndxi2dxi2 + aInvJac(:,3).*d2Ndxi3dxi2 ) + ...
              aInvJac(:,3).*( aInvJac(:,1).*d2Ndxi1dxi3 + aInvJac(:,2).*d2Ndxi2dxi3 + aInvJac(:,3).*d2Ndxi3dxi3 );
    elseif( i_eval==33 )
      vBase = aInvJac(:,4).*( aInvJac(:,4).*d2Ndxi1dxi1 + aInvJac(:,5).*d2Ndxi2dxi1 + aInvJac(:,6).*d2Ndxi3dxi1 ) + ...
              aInvJac(:,5).*( aInvJac(:,4).*d2Ndxi1dxi2 + aInvJac(:,5).*d2Ndxi2dxi2 + aInvJac(:,6).*d2Ndxi3dxi2 ) + ...
              aInvJac(:,6).*( aInvJac(:,4).*d2Ndxi1dxi3 + aInvJac(:,5).*d2Ndxi2dxi3 + aInvJac(:,6).*d2Ndxi3dxi3 );
    else
      vBase = aInvJac(:,4).*( aInvJac(:,1).*d2Ndxi1dxi1 + aInvJac(:,2).*d2Ndxi2dxi1 + aInvJac(:,3).*d2Ndxi3dxi1 ) + ...
              aInvJac(:,5).*( aInvJac(:,1).*d2Ndxi1dxi2 + aInvJac(:,2).*d2Ndxi2dxi2 + aInvJac(:,3).*d2Ndxi3dxi2 ) + ...
              aInvJac(:,6).*( aInvJac(:,1).*d2Ndxi1dxi3 + aInvJac(:,2).*d2Ndxi2dxi3 + aInvJac(:,3).*d2Ndxi3dxi3 );
    end

  otherwise
    vBase = 0;

end
