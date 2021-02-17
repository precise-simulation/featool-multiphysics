function [ vBase, nLDof, xLDof, sfun ] = sf_tet_P5( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_TET_P5 Fifth order Lagrange shape functions for tetrahedra (P5).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_TET_P5( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fifth order P5 Lagrange shape functions on 3D tetrahedral elements
%   with values defined in the nodes, edges, faces, and cell center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 3              Number of space dimensions
%       n_vert      scalar: 4              Number of vertices per cell
%       i_dof       scalar: 1-56           Local basis function to evaluate
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


nLDof = [4 24 24 4];
xLDof = [1 0 0 0 4/5   0 1/5 4/5   0   0 3/5   0 2/5 3/5   0   0 2/5   0 3/5 2/5   0   0 1/5   0 4/5 1/5   0   0 3/5 3/5   0 1/5 2/5 2/5   0 2/5 1/5 1/5   0 3/5 1/5 1/5   0 2/5 1/5 1/5   0 1/5 2/5 2/5   0 1/5 2/5 1/5 1/5 1/5;
         0 1 0 0 1/5 4/5   0   0 4/5   0 2/5 3/5   0   0 3/5   0 3/5 2/5   0   0 2/5   0 4/5 1/5   0   0 1/5   0 1/5 1/5 3/5   0 2/5 2/5 2/5   0 3/5 3/5 1/5   0 2/5 2/5 1/5   0 1/5 1/5 1/5   0 1/5 1/5 2/5   0 1/5 2/5 1/5 1/5;
         0 0 1 0   0 1/5 4/5   0   0 4/5   0 2/5 3/5   0   0 3/5   0 3/5 2/5   0   0 2/5   0 4/5 1/5   0   0 1/5 1/5   0 1/5 3/5 1/5   0 2/5 2/5 1/5   0 3/5 1/5 2/5   0 2/5 1/5 3/5   0 1/5 1/5 2/5   0 1/5 2/5 1/5 1/5 2/5 1/5;
         0 0 0 1   0   0   0 1/5 1/5 1/5   0   0   0 2/5 2/5 2/5   0   0   0 3/5 3/5 3/5   0   0   0 4/5 4/5 4/5   0 1/5 1/5 1/5   0 1/5 1/5 1/5   0 1/5 1/5 1/5   0 2/5 2/5 2/5   0 3/5 3/5 3/5   0 2/5 2/5 2/5 1/5 1/5 1/5 2/5];
sfun  = 'sf_tet_P5';

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
        vBase = (xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3)*(5*xi(4) - 4))/24;
      case 5
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 6
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 7
        vBase = (25*xi(1)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 8
        vBase = (25*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 9
        vBase = (25*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 10
        vBase = (25*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 11
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(2) - 1))/12;
      case 12
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(3) - 1))/12;
      case 13
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
      case 14
        vBase = (25*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(4) - 1))/12;
      case 15
        vBase = (25*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(4) - 1))/12;
      case 16
        vBase = (25*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(4) - 1))/12;
      case 17
        vBase = (25*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1)*(5*xi(2) - 2))/12;
      case 18
        vBase = (25*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
      case 19
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(3) - 1))/12;
      case 20
        vBase = (25*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
      case 21
        vBase = (25*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
      case 22
        vBase = (25*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
      case 23
        vBase = (25*xi(1)*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 24
        vBase = (25*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 25
        vBase = (25*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 26
        vBase = (25*xi(1)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
      case 27
        vBase = (25*xi(2)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
      case 28
        vBase = (25*xi(3)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
      case 29
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 30
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 31
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 32
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 33
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
      case 34
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
      case 35
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
      case 36
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
      case 37
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 38
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 39
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 40
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 41
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
      case 42
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1))/4;
      case 43
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1))/4;
      case 44
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1))/4;
      case 45
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 46
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
      case 47
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
      case 48
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
      case 49
        vBase = (125*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
      case 50
        vBase = (125*xi(1)*xi(2)*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1))/4;
      case 51
        vBase = (125*xi(2)*xi(3)*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1))/4;
      case 52
        vBase = (125*xi(1)*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1))/4;
      case 53
        vBase = (625*xi(1)*xi(2)*xi(3)*xi(4)*(5*xi(1) - 1))/2;
      case 54
        vBase = (625*xi(1)*xi(2)*xi(3)*xi(4)*(5*xi(2) - 1))/2;
      case 55
        vBase = (625*xi(1)*xi(2)*xi(3)*xi(4)*(5*xi(3) - 1))/2;
      case 56
        vBase = (625*xi(1)*xi(2)*xi(3)*xi(4)*(5*xi(4) - 1))/2;
    end

  case {2,3,4}

    switch i_dof

      case 1
        dNdxi1 = (3125*xi(1)^4)/24 - (625*xi(1)^3)/3 + (875*xi(1)^2)/8 - (125*xi(1))/6 + 1;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 2
        dNdxi1 = 0;
        dNdxi2 = (3125*xi(2)^4)/24 - (625*xi(2)^3)/3 + (875*xi(2)^2)/8 - (125*xi(2))/6 + 1;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 3
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (3125*xi(3)^4)/24 - (625*xi(3)^3)/3 + (875*xi(3)^2)/8 - (125*xi(3))/6 + 1;
        dNdxi4 = 0;
      case 4
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (3125*xi(4)^4)/24 - (625*xi(4)^3)/3 + (875*xi(4)^2)/8 - (125*xi(4))/6 + 1;
      case 5
        dNdxi1 = (25*xi(2)*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 6
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
        dNdxi4 = 0;
      case 7
        dNdxi1 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        dNdxi4 = 0;
      case 8
        dNdxi1 = (25*xi(4)*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
      case 9
        dNdxi1 = 0;
        dNdxi2 = (25*xi(4)*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
      case 10
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(4)*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        dNdxi4 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
      case 11
        dNdxi1 = (25*xi(2)*(5*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(10*xi(2) - 1))/12;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 12
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(5*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(10*xi(3) - 1))/12;
        dNdxi4 = 0;
      case 13
        dNdxi1 = (25*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        dNdxi4 = 0;
      case 14
        dNdxi1 = (25*xi(4)*(5*xi(4) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(10*xi(4) - 1))/12;
      case 15
        dNdxi1 = 0;
        dNdxi2 = (25*xi(4)*(5*xi(4) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(10*xi(4) - 1))/12;
      case 16
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(4)*(5*xi(4) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        dNdxi4 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(10*xi(4) - 1))/12;
      case 17
        dNdxi1 = (25*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1)*(5*xi(2) - 2))/12;
        dNdxi2 = (25*xi(1)*(5*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 18
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1)*(5*xi(3) - 2))/12;
        dNdxi3 = (25*xi(2)*(5*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        dNdxi4 = 0;
      case 19
        dNdxi1 = (25*xi(3)*(5*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(10*xi(3) - 1))/12;
        dNdxi4 = 0;
      case 20
        dNdxi1 = (25*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(1)*(5*xi(1) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
      case 21
        dNdxi1 = 0;
        dNdxi2 = (25*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(2)*(5*xi(2) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
      case 22
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1)*(5*xi(4) - 2))/12;
        dNdxi4 = (25*xi(3)*(5*xi(3) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
      case 23
        dNdxi1 = (25*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2)*(5*xi(2) - 3))/24;
        dNdxi2 = (25*xi(1)*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        dNdxi3 = 0;
        dNdxi4 = 0;
      case 24
        dNdxi1 = 0;
        dNdxi2 = (25*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2)*(5*xi(3) - 3))/24;
        dNdxi3 = (25*xi(2)*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        dNdxi4 = 0;
      case 25
        dNdxi1 = (25*xi(3)*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2)*(5*xi(1) - 3))/24;
        dNdxi4 = 0;
      case 26
        dNdxi1 = (25*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
        dNdxi2 = 0;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(1)*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
      case 27
        dNdxi1 = 0;
        dNdxi2 = (25*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
        dNdxi3 = 0;
        dNdxi4 = (25*xi(2)*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
      case 28
        dNdxi1 = 0;
        dNdxi2 = 0;
        dNdxi3 = (25*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2)*(5*xi(4) - 3))/24;
        dNdxi4 = (25*xi(3)*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
      case 29
        dNdxi1 = (125*xi(2)*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        dNdxi4 = 0;
      case 30
        dNdxi1 = (125*xi(2)*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        dNdxi2 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 31
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        dNdxi3 = (125*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        dNdxi4 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 32
        dNdxi1 = (125*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        dNdxi4 = (125*xi(1)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 33
        dNdxi1 = (125*xi(2)*xi(3)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
        dNdxi4 = 0;
      case 34
        dNdxi1 = (125*xi(2)*xi(4)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1))/4;
      case 35
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        dNdxi3 = (125*xi(2)*xi(4)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        dNdxi4 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
      case 36
        dNdxi1 = (125*xi(3)*xi(4)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        dNdxi4 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
      case 37
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        dNdxi4 = 0;
      case 38
        dNdxi1 = (125*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
      case 39
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi3 = (125*xi(2)*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        dNdxi4 = (125*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
      case 40
        dNdxi1 = (125*xi(3)*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        dNdxi4 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
      case 41
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        dNdxi4 = 0;
      case 42
        dNdxi1 = (125*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
      case 43
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1))/4;
        dNdxi3 = (125*xi(2)*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        dNdxi4 = (125*xi(2)*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
      case 44
        dNdxi1 = (125*xi(3)*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1))/4;
        dNdxi4 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
      case 45
        dNdxi1 = (125*xi(2)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        dNdxi3 = (125*xi(1)*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        dNdxi4 = 0;
      case 46
        dNdxi1 = (125*xi(2)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi2 = (125*xi(1)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
      case 47
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi3 = (125*xi(2)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi4 = (125*xi(2)*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
      case 48
        dNdxi1 = (125*xi(3)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        dNdxi4 = (125*xi(1)*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
      case 49
        dNdxi1 = (125*xi(2)*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/4;
        dNdxi3 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        dNdxi4 = 0;
      case 50
        dNdxi1 = (125*xi(2)*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        dNdxi2 = (125*xi(1)*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1))/4;
        dNdxi3 = 0;
        dNdxi4 = (125*xi(1)*xi(2)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
      case 51
        dNdxi1 = 0;
        dNdxi2 = (125*xi(3)*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        dNdxi3 = (125*xi(2)*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1))/4;
        dNdxi4 = (125*xi(2)*xi(3)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
      case 52
        dNdxi1 = (125*xi(3)*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1))/4;
        dNdxi2 = 0;
        dNdxi3 = (125*xi(1)*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        dNdxi4 = (125*xi(1)*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
      case 53
        dNdxi1 = (625*xi(2)*xi(3)*xi(4)*(10*xi(1) - 1))/2;
        dNdxi2 = (625*xi(1)*xi(3)*xi(4)*(5*xi(1) - 1))/2;
        dNdxi3 = (625*xi(1)*xi(2)*xi(4)*(5*xi(1) - 1))/2;
        dNdxi4 = (625*xi(1)*xi(2)*xi(3)*(5*xi(1) - 1))/2;
      case 54
        dNdxi1 = (625*xi(2)*xi(3)*xi(4)*(5*xi(2) - 1))/2;
        dNdxi2 = (625*xi(1)*xi(3)*xi(4)*(10*xi(2) - 1))/2;
        dNdxi3 = (625*xi(1)*xi(2)*xi(4)*(5*xi(2) - 1))/2;
        dNdxi4 = (625*xi(1)*xi(2)*xi(3)*(5*xi(2) - 1))/2;
      case 55
        dNdxi1 = (625*xi(2)*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        dNdxi2 = (625*xi(1)*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        dNdxi3 = (625*xi(1)*xi(2)*xi(4)*(10*xi(3) - 1))/2;
        dNdxi4 = (625*xi(1)*xi(2)*xi(3)*(5*xi(3) - 1))/2;
      case 56
        dNdxi1 = (625*xi(2)*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        dNdxi2 = (625*xi(1)*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        dNdxi3 = (625*xi(1)*xi(2)*xi(4)*(5*xi(4) - 1))/2;
        dNdxi4 = (625*xi(1)*xi(2)*xi(3)*(10*xi(4) - 1))/2;
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
        d2Ndxi1dxi1 = (125*(5*xi(1) - 2)*(10*xi(1)^2 - 8*xi(1) + 1))/12;
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
        d2Ndxi2dxi2 = (125*(5*xi(2) - 2)*(10*xi(2)^2 - 8*xi(2) + 1))/12;
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
        d2Ndxi3dxi3 = (125*(5*xi(3) - 2)*(10*xi(3)^2 - 8*xi(3) + 1))/12;
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
        d2Ndxi4dxi4 = (125*(5*xi(4) - 2)*(10*xi(4)^2 - 8*xi(4) + 1))/12;

      case 5
        d2Ndxi1dxi1 = (125*xi(2)*(150*xi(1)^2 - 90*xi(1) + 11))/12;
        d2Ndxi2dxi1 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
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
        d2Ndxi2dxi2 = (125*xi(3)*(150*xi(2)^2 - 90*xi(2) + 11))/12;
        d2Ndxi3dxi2 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 7
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(150*xi(3)^2 - 90*xi(3) + 11))/12;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 8
        d2Ndxi1dxi1 = (125*xi(4)*(150*xi(1)^2 - 90*xi(1) + 11))/12;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 9
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(4)*(150*xi(2)^2 - 90*xi(2) + 11))/12;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
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
        d2Ndxi3dxi3 = (125*xi(4)*(150*xi(3)^2 - 90*xi(3) + 11))/12;
        d2Ndxi4dxi3 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi4dxi4 = 0;

      case 11
        d2Ndxi1dxi1 = (125*xi(2)*(5*xi(1) - 1)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi1 = (25*(10*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(2) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
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
        d2Ndxi2dxi2 = (125*xi(3)*(5*xi(2) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi2 = (25*(10*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(3) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 13
        d2Ndxi1dxi1 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(1) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 14
        d2Ndxi1dxi1 = (125*xi(4)*(5*xi(1) - 1)*(5*xi(4) - 1))/2;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (25*(10*xi(4) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (25*(10*xi(4) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;

      case 15
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(4)*(5*xi(2) - 1)*(5*xi(4) - 1))/2;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (25*(10*xi(4) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (25*(10*xi(4) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;

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
        d2Ndxi3dxi3 = (125*xi(4)*(5*xi(3) - 1)*(5*xi(4) - 1))/2;
        d2Ndxi4dxi3 = (25*(10*xi(4) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (25*(10*xi(4) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi4dxi4 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;

      case 17
        d2Ndxi1dxi1 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi2dxi1 = (25*(10*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(1) - 1)*(75*xi(2)^2 - 30*xi(2) + 2))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(2) - 1))/2;
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
        d2Ndxi2dxi2 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi3dxi2 = (25*(10*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(2) - 1)*(75*xi(3)^2 - 30*xi(3) + 2))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 19
        d2Ndxi1dxi1 = (125*xi(3)*(5*xi(1) - 1)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(3) - 1)*(75*xi(1)^2 - 30*xi(1) + 2))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 20
        d2Ndxi1dxi1 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (25*(10*xi(1) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (25*(10*xi(1) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(4) - 1))/2;

      case 21
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (25*(10*xi(2) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (25*(10*xi(2) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(4) - 1))/2;

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
        d2Ndxi3dxi3 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi4dxi3 = (25*(10*xi(3) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (25*(10*xi(3) - 1)*(75*xi(4)^2 - 30*xi(4) + 2))/12;
        d2Ndxi4dxi4 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(4) - 1))/2;

      case 23
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (25*(10*xi(2) - 3)*(25*xi(2)^2 - 15*xi(2) + 1))/12;
        d2Ndxi2dxi2 = (125*xi(1)*(150*xi(2)^2 - 90*xi(2) + 11))/12;
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

      case 24
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (25*(10*xi(3) - 3)*(25*xi(3)^2 - 15*xi(3) + 1))/12;
        d2Ndxi3dxi3 = (125*xi(2)*(150*xi(3)^2 - 90*xi(3) + 11))/12;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 25
        d2Ndxi1dxi1 = (125*xi(3)*(150*xi(1)^2 - 90*xi(1) + 11))/12;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (25*(10*xi(1) - 3)*(25*xi(1)^2 - 15*xi(1) + 1))/12;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 26
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(1)*(150*xi(4)^2 - 90*xi(4) + 11))/12;

      case 27
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (125*xi(2)*(150*xi(4)^2 - 90*xi(4) + 11))/12;

      case 28
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
        d2Ndxi4dxi3 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (25*(10*xi(4) - 3)*(25*xi(4)^2 - 15*xi(4) + 1))/12;
        d2Ndxi4dxi4 = (125*xi(3)*(150*xi(4)^2 - 90*xi(4) + 11))/12;

      case 29
        d2Ndxi1dxi1 = 625*xi(2)*xi(3)*(5*xi(1) - 1);
        d2Ndxi2dxi1 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 30
        d2Ndxi1dxi1 = 625*xi(2)*xi(4)*(5*xi(1) - 1);
        d2Ndxi2dxi1 = (125*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi1dxi2 = (125*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi4 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 31
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 625*xi(3)*xi(4)*(5*xi(2) - 1);
        d2Ndxi3dxi2 = (125*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi4dxi2 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi4 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi4dxi4 = 0;

      case 32
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi1 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 625*xi(1)*xi(4)*(5*xi(3) - 1);
        d2Ndxi4dxi3 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi1dxi4 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi4 = 0;

      case 33
        d2Ndxi1dxi1 = (625*xi(2)*xi(3)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(3)*(5*xi(1) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 34
        d2Ndxi1dxi1 = (625*xi(2)*xi(4)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(4)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(4)*(10*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(4)*(5*xi(1) - 1))/2;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(10*xi(1) - 1)*(5*xi(2) - 1))/4;
        d2Ndxi2dxi4 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(2) - 1))/4;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 35
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (625*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(4)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi2 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(2)*xi(4)*(5*xi(2) - 1))/2;
        d2Ndxi4dxi3 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi3dxi4 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi4 = 0;

      case 36
        d2Ndxi1dxi1 = (625*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi1 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (625*xi(1)*xi(4)*(5*xi(1) - 1))/2;
        d2Ndxi4dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi1dxi4 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi4 = 0;

      case 37
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi2dxi2 = 625*xi(1)*xi(3)*(5*xi(2) - 1);
        d2Ndxi3dxi2 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 38
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi1dxi2 = (125*xi(4)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi2dxi2 = 625*xi(1)*xi(4)*(5*xi(2) - 1);
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(5*xi(2) - 1)*(5*xi(2) - 2))/6;
        d2Ndxi2dxi4 = (125*xi(1)*(75*xi(2)^2 - 30*xi(2) + 2))/6;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 39
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi2 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi3dxi3 = 625*xi(2)*xi(4)*(5*xi(3) - 1);
        d2Ndxi4dxi3 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi3dxi4 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi4 = 0;

      case 40
        d2Ndxi1dxi1 = 625*xi(3)*xi(4)*(5*xi(1) - 1);
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi4dxi1 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi1dxi4 = (125*xi(3)*(75*xi(1)^2 - 30*xi(1) + 2))/6;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(5*xi(1) - 1)*(5*xi(1) - 2))/6;
        d2Ndxi4dxi4 = 0;

      case 41
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(2) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(1)*xi(2)*(5*xi(2) - 1))/2;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 42
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi2dxi2 = (625*xi(1)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi2dxi4 = (125*xi(1)*(10*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (625*xi(1)*xi(2)*(5*xi(2) - 1))/2;

      case 43
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi4dxi2 = (125*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(2)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi4dxi3 = (125*xi(2)*(10*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi3dxi4 = (125*xi(2)*(10*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi4dxi4 = (625*xi(2)*xi(3)*(5*xi(3) - 1))/2;

      case 44
        d2Ndxi1dxi1 = (625*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi4dxi1 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi4 = (125*xi(3)*(10*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi4dxi4 = (625*xi(1)*xi(3)*(5*xi(1) - 1))/2;

      case 45
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi3dxi1 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(5*xi(3) - 1)*(5*xi(3) - 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi2dxi3 = (125*xi(1)*(75*xi(3)^2 - 30*xi(3) + 2))/6;
        d2Ndxi3dxi3 = 625*xi(1)*xi(2)*(5*xi(3) - 1);
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 46
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi2 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi2dxi4 = (125*xi(1)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 625*xi(1)*xi(2)*(5*xi(4) - 1);

      case 47
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi4dxi2 = (125*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(2)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi3dxi4 = (125*xi(2)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi4dxi4 = 625*xi(2)*xi(3)*(5*xi(4) - 1);

      case 48
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi4dxi1 = (125*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(5*xi(4) - 1)*(5*xi(4) - 2))/6;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(1)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi1dxi4 = (125*xi(3)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(75*xi(4)^2 - 30*xi(4) + 2))/6;
        d2Ndxi4dxi4 = 625*xi(1)*xi(3)*(5*xi(4) - 1);

      case 49
        d2Ndxi1dxi1 = (625*xi(2)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi3dxi1 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = (125*xi(3)*(10*xi(1) - 1)*(5*xi(3) - 1))/4;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi2dxi3 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(3) - 1))/4;
        d2Ndxi3dxi3 = (625*xi(1)*xi(2)*(5*xi(1) - 1))/2;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = 0;

      case 50
        d2Ndxi1dxi1 = (625*xi(2)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi2dxi1 = (125*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi2 = (125*xi(4)*(10*xi(1) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = 0;
        d2Ndxi1dxi4 = (125*xi(2)*(10*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi2dxi4 = (125*xi(1)*(5*xi(1) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi3dxi4 = 0;
        d2Ndxi4dxi4 = (625*xi(1)*xi(2)*(5*xi(1) - 1))/2;

      case 51
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = 0;
        d2Ndxi4dxi1 = 0;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = (625*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi3dxi2 = (125*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi4dxi2 = (125*xi(3)*(10*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi3 = 0;
        d2Ndxi2dxi3 = (125*xi(4)*(10*xi(2) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi4 = 0;
        d2Ndxi2dxi4 = (125*xi(3)*(10*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi3dxi4 = (125*xi(2)*(5*xi(2) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi4dxi4 = (625*xi(2)*xi(3)*(5*xi(2) - 1))/2;

      case 52
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi3dxi1 = (125*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi4dxi1 = (125*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 0;
        d2Ndxi4dxi2 = 0;
        d2Ndxi1dxi3 = (125*xi(4)*(10*xi(3) - 1)*(5*xi(4) - 1))/4;
        d2Ndxi2dxi3 = 0;
        d2Ndxi3dxi3 = (625*xi(1)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi4dxi3 = (125*xi(1)*(10*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi1dxi4 = (125*xi(3)*(5*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi2dxi4 = 0;
        d2Ndxi3dxi4 = (125*xi(1)*(10*xi(3) - 1)*(10*xi(4) - 1))/4;
        d2Ndxi4dxi4 = (625*xi(1)*xi(3)*(5*xi(3) - 1))/2;

      case 53
        d2Ndxi1dxi1 = 3125*xi(2)*xi(3)*xi(4);
        d2Ndxi2dxi1 = (625*xi(3)*xi(4)*(10*xi(1) - 1))/2;
        d2Ndxi3dxi1 = (625*xi(2)*xi(4)*(10*xi(1) - 1))/2;
        d2Ndxi4dxi1 = (625*xi(2)*xi(3)*(10*xi(1) - 1))/2;
        d2Ndxi1dxi2 = (625*xi(3)*xi(4)*(10*xi(1) - 1))/2;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (625*xi(1)*xi(4)*(5*xi(1) - 1))/2;
        d2Ndxi4dxi2 = (625*xi(1)*xi(3)*(5*xi(1) - 1))/2;
        d2Ndxi1dxi3 = (625*xi(2)*xi(4)*(10*xi(1) - 1))/2;
        d2Ndxi2dxi3 = (625*xi(1)*xi(4)*(5*xi(1) - 1))/2;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (625*xi(1)*xi(2)*(5*xi(1) - 1))/2;
        d2Ndxi1dxi4 = (625*xi(2)*xi(3)*(10*xi(1) - 1))/2;
        d2Ndxi2dxi4 = (625*xi(1)*xi(3)*(5*xi(1) - 1))/2;
        d2Ndxi3dxi4 = (625*xi(1)*xi(2)*(5*xi(1) - 1))/2;
        d2Ndxi4dxi4 = 0;

      case 54
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (625*xi(3)*xi(4)*(10*xi(2) - 1))/2;
        d2Ndxi3dxi1 = (625*xi(2)*xi(4)*(5*xi(2) - 1))/2;
        d2Ndxi4dxi1 = (625*xi(2)*xi(3)*(5*xi(2) - 1))/2;
        d2Ndxi1dxi2 = (625*xi(3)*xi(4)*(10*xi(2) - 1))/2;
        d2Ndxi2dxi2 = 3125*xi(1)*xi(3)*xi(4);
        d2Ndxi3dxi2 = (625*xi(1)*xi(4)*(10*xi(2) - 1))/2;
        d2Ndxi4dxi2 = (625*xi(1)*xi(3)*(10*xi(2) - 1))/2;
        d2Ndxi1dxi3 = (625*xi(2)*xi(4)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi3 = (625*xi(1)*xi(4)*(10*xi(2) - 1))/2;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (625*xi(1)*xi(2)*(5*xi(2) - 1))/2;
        d2Ndxi1dxi4 = (625*xi(2)*xi(3)*(5*xi(2) - 1))/2;
        d2Ndxi2dxi4 = (625*xi(1)*xi(3)*(10*xi(2) - 1))/2;
        d2Ndxi3dxi4 = (625*xi(1)*xi(2)*(5*xi(2) - 1))/2;
        d2Ndxi4dxi4 = 0;

      case 55
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (625*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi1 = (625*xi(2)*xi(4)*(10*xi(3) - 1))/2;
        d2Ndxi4dxi1 = (625*xi(2)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi1dxi2 = (625*xi(3)*xi(4)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (625*xi(1)*xi(4)*(10*xi(3) - 1))/2;
        d2Ndxi4dxi2 = (625*xi(1)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi1dxi3 = (625*xi(2)*xi(4)*(10*xi(3) - 1))/2;
        d2Ndxi2dxi3 = (625*xi(1)*xi(4)*(10*xi(3) - 1))/2;
        d2Ndxi3dxi3 = 3125*xi(1)*xi(2)*xi(4);
        d2Ndxi4dxi3 = (625*xi(1)*xi(2)*(10*xi(3) - 1))/2;
        d2Ndxi1dxi4 = (625*xi(2)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi2dxi4 = (625*xi(1)*xi(3)*(5*xi(3) - 1))/2;
        d2Ndxi3dxi4 = (625*xi(1)*xi(2)*(10*xi(3) - 1))/2;
        d2Ndxi4dxi4 = 0;

      case 56
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = (625*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi3dxi1 = (625*xi(2)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi4dxi1 = (625*xi(2)*xi(3)*(10*xi(4) - 1))/2;
        d2Ndxi1dxi2 = (625*xi(3)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = (625*xi(1)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi4dxi2 = (625*xi(1)*xi(3)*(10*xi(4) - 1))/2;
        d2Ndxi1dxi3 = (625*xi(2)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi2dxi3 = (625*xi(1)*xi(4)*(5*xi(4) - 1))/2;
        d2Ndxi3dxi3 = 0;
        d2Ndxi4dxi3 = (625*xi(1)*xi(2)*(10*xi(4) - 1))/2;
        d2Ndxi1dxi4 = (625*xi(2)*xi(3)*(10*xi(4) - 1))/2;
        d2Ndxi2dxi4 = (625*xi(1)*xi(3)*(10*xi(4) - 1))/2;
        d2Ndxi3dxi4 = (625*xi(1)*xi(2)*(10*xi(4) - 1))/2;
        d2Ndxi4dxi4 = 3125*xi(1)*xi(2)*xi(3);

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
