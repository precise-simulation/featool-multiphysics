function [ vBase, nLDof, xLDof, sfun ] = sf_hex_Q2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_HEX_Q2 Triquadratic conforming shape function for hexahedrons (Q2).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_HEX_Q2( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming triquadratic Q2 shape functions on hexahedrons with values defined
%   in the nodes, edges, faces, and cell center. XI is [-1..1]^3 reference coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  3             Number of space dimensions
%       n_vert      scalar:  8             Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim]               Local coordinates of evaluation point
%       aInvJac     [n,n_sdim*n_sdim]      Inverse of transformation Jacobian
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
%   See also SFLAG2, SF_HEX_Q1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [8 12 6 1];
xLDof = [-1  1  1 -1  -1  1  1 -1   0  1  0 -1  -1  1  1 -1   0  1  0 -1   0  0  1 0 -1 0  0; ...
         -1 -1  1  1  -1 -1  1  1  -1  0  1  0  -1 -1  1  1  -1  0  1  0   0 -1  0 1  0 0  0; ...
         -1 -1 -1 -1   1  1  1  1  -1 -1 -1 -1   0  0  0  0   1  1  1  1  -1  0  0 0  0 1  0];
sfun  = 'sf_hex_Q2';


switch i_eval   % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = -xi(1)*(1-xi(1))*xi(2)*(1-xi(2))*xi(3)*(1-xi(3))/8;
      case 2
        vBase =  xi(1)*(1+xi(1))*xi(2)*(1-xi(2))*xi(3)*(1-xi(3))/8;
      case 3
        vBase = -xi(1)*(1+xi(1))*xi(2)*(1+xi(2))*xi(3)*(1-xi(3))/8;
      case 4
        vBase =  xi(1)*(1-xi(1))*xi(2)*(1+xi(2))*xi(3)*(1-xi(3))/8;
      case 5
        vBase =  xi(1)*(1-xi(1))*xi(2)*(1-xi(2))*xi(3)*(1+xi(3))/8;
      case 6
        vBase = -xi(1)*(1+xi(1))*xi(2)*(1-xi(2))*xi(3)*(1+xi(3))/8;
      case 7
        vBase =  xi(1)*(1+xi(1))*xi(2)*(1+xi(2))*xi(3)*(1+xi(3))/8;
      case 8
        vBase = -xi(1)*(1-xi(1))*xi(2)*(1+xi(2))*xi(3)*(1+xi(3))/8;
      case 9
        vBase =  (2-2*xi(1)^2)*xi(2)*(1-xi(2))*xi(3)*(1-xi(3))/8;
      case 10
        vBase = -xi(1)*(1+xi(1))*(2-2*xi(2)^2)*xi(3)*(1-xi(3))/8;
      case 11
        vBase = -(2-2*xi(1)^2)*xi(2)*(1+xi(2))*xi(3)*(1-xi(3))/8;
      case 12
        vBase =  xi(1)*(1-xi(1))*(2-2*xi(2)^2)*xi(3)*(1-xi(3))/8;
      case 13
        vBase =  xi(1)*(1-xi(1))*xi(2)*(1-xi(2))*(2-2*xi(3)^2)/8;
      case 14
        vBase = -xi(1)*(1+xi(1))*xi(2)*(1-xi(2))*(2-2*xi(3)^2)/8;
      case 15
        vBase =  xi(1)*(1+xi(1))*xi(2)*(1+xi(2))*(2-2*xi(3)^2)/8;
      case 16
        vBase = -xi(1)*(1-xi(1))*xi(2)*(1+xi(2))*(2-2*xi(3)^2)/8;
      case 17
        vBase = -(2-2*xi(1)^2)*xi(2)*(1-xi(2))*xi(3)*(1+xi(3))/8;
      case 18
        vBase =  xi(1)*(1+xi(1))*(2-2*xi(2)^2)*xi(3)*(1+xi(3))/8;
      case 19
        vBase =  (2-2*xi(1)^2)*xi(2)*(1+xi(2))*xi(3)*(1+xi(3))/8;
      case 20
        vBase = -xi(1)*(1-xi(1))*(2-2*xi(2)^2)*xi(3)*(1+xi(3))/8;
      case 21
        vBase = -(2-2*xi(1)^2)*(2-2*xi(2)^2)*xi(3)*(1-xi(3))/8;
      case 22
        vBase = -(2-2*xi(1)^2)*xi(2)*(1-xi(2))*(2-2*xi(3)^2)/8;
      case 23
        vBase =  xi(1)*(1+xi(1))*(2-2*xi(2)^2)*(2-2*xi(3)^2)/8;
      case 24
        vBase =  (2-2*xi(1)^2)*xi(2)*(1+xi(2))*(2-2*xi(3)^2)/8;
      case 25
        vBase = -xi(1)*(1-xi(1))*(2-2*xi(2)^2)*(2-2*xi(3)^2)/8;
      case 26
        vBase =  (2-2*xi(1)^2)*(2-2*xi(2)^2)*xi(3)*(1+xi(3))/8;
      case 27
        vBase =  (2-2*xi(1)^2)*(2-2*xi(2)^2)*(2-2*xi(3)^2)/8;

    end

  case {2,3,4}   % Evaluation of first order derivatives.

    switch i_dof   % Basis function to evaluate.
      case 1
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)-1)*(xi(2)-1)*(xi(3)-1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)-1)*(xi(1)-1)*(xi(3)-1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)-1)*(xi(1)-1)*(xi(2)-1))/8;
      case 2
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)+1)*(xi(2)-1)*(xi(3)-1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)-1)*(xi(1)+1)*(xi(3)-1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)-1)*(xi(1)+1)*(xi(2)-1))/8;
      case 3
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)+1)*(xi(2)+1)*(xi(3)-1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)+1)*(xi(1)+1)*(xi(3)-1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)-1)*(xi(1)+1)*(xi(2)+1))/8;
      case 4
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)-1)*(xi(2)+1)*(xi(3)-1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)+1)*(xi(1)-1)*(xi(3)-1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)-1)*(xi(1)-1)*(xi(2)+1))/8;
      case 5
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)-1)*(xi(2)-1)*(xi(3)+1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)-1)*(xi(1)-1)*(xi(3)+1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)+1)*(xi(1)-1)*(xi(2)-1))/8;
      case 6
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)+1)*(xi(2)-1)*(xi(3)+1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)-1)*(xi(1)+1)*(xi(3)+1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)+1)*(xi(1)+1)*(xi(2)-1))/8;
      case 7
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)+1)*(xi(2)+1)*(xi(3)+1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)+1)*(xi(1)+1)*(xi(3)+1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)+1)*(xi(1)+1)*(xi(2)+1))/8;
      case 8
        dNdxi1 = (xi(2)*xi(3)*(2*xi(1)-1)*(xi(2)+1)*(xi(3)+1))/8;
        dNdxi2 = (xi(1)*xi(3)*(2*xi(2)+1)*(xi(1)-1)*(xi(3)+1))/8;
        dNdxi3 = (xi(1)*xi(2)*(2*xi(3)+1)*(xi(1)-1)*(xi(2)+1))/8;
      case 9
        dNdxi1 = -(xi(1)*xi(2)*xi(3)*(xi(2)-1)*(xi(3)-1))/2;
        dNdxi2 = -(xi(3)*(2*xi(2)-1)*(xi(1)^2-1)*(xi(3)-1))/4;
        dNdxi3 = -(xi(2)*(2*xi(3)-1)*(xi(1)^2-1)*(xi(2)-1))/4;
      case 10
        dNdxi1 = -(xi(3)*(2*xi(1)+1)*(xi(2)^2-1)*(xi(3)-1))/4;
        dNdxi2 = -(xi(1)*xi(2)*xi(3)*(xi(1)+1)*(xi(3)-1))/2;
        dNdxi3 = -(xi(1)*(2*xi(3)-1)*(xi(2)^2-1)*(xi(1)+1))/4;
      case 11
        dNdxi1 = -(xi(1)*xi(2)*xi(3)*(xi(2)+1)*(xi(3)-1))/2;
        dNdxi2 = -(xi(3)*(2*xi(2)+1)*(xi(1)^2-1)*(xi(3)-1))/4;
        dNdxi3 = -(xi(2)*(2*xi(3)-1)*(xi(1)^2-1)*(xi(2)+1))/4;
      case 12
        dNdxi1 = -(xi(3)*(2*xi(1)-1)*(xi(2)^2-1)*(xi(3)-1))/4;
        dNdxi2 = -(xi(1)*xi(2)*xi(3)*(xi(1)-1)*(xi(3)-1))/2;
        dNdxi3 = -(xi(1)*(2*xi(3)-1)*(xi(2)^2-1)*(xi(1)-1))/4;
      case 13
        dNdxi1 = -(xi(2)*(2*xi(1)-1)*(xi(3)^2-1)*(xi(2)-1))/4;
        dNdxi2 = -(xi(1)*(2*xi(2)-1)*(xi(3)^2-1)*(xi(1)-1))/4;
        dNdxi3 = -(xi(1)*xi(2)*xi(3)*(xi(1)-1)*(xi(2)-1))/2;
      case 14
        dNdxi1 = -(xi(2)*(2*xi(1)+1)*(xi(3)^2-1)*(xi(2)-1))/4;
        dNdxi2 = -(xi(1)*(2*xi(2)-1)*(xi(3)^2-1)*(xi(1)+1))/4;
        dNdxi3 = -(xi(1)*xi(2)*xi(3)*(xi(1)+1)*(xi(2)-1))/2;
      case 15
        dNdxi1 = -(xi(2)*(2*xi(1)+1)*(xi(3)^2-1)*(xi(2)+1))/4;
        dNdxi2 = -(xi(1)*(2*xi(2)+1)*(xi(3)^2-1)*(xi(1)+1))/4;
        dNdxi3 = -(xi(1)*xi(2)*xi(3)*(xi(1)+1)*(xi(2)+1))/2;
      case 16
        dNdxi1 = -(xi(2)*(2*xi(1)-1)*(xi(3)^2-1)*(xi(2)+1))/4;
        dNdxi2 = -(xi(1)*(2*xi(2)+1)*(xi(3)^2-1)*(xi(1)-1))/4;
        dNdxi3 = -(xi(1)*xi(2)*xi(3)*(xi(1)-1)*(xi(2)+1))/2;
      case 17
        dNdxi1 = -(xi(1)*xi(2)*xi(3)*(xi(2)-1)*(xi(3)+1))/2;
        dNdxi2 = -(xi(3)*(2*xi(2)-1)*(xi(1)^2-1)*(xi(3)+1))/4;
        dNdxi3 = -(xi(2)*(2*xi(3)+1)*(xi(1)^2-1)*(xi(2)-1))/4;
      case 18
        dNdxi1 = -(xi(3)*(2*xi(1)+1)*(xi(2)^2-1)*(xi(3)+1))/4;
        dNdxi2 = -(xi(1)*xi(2)*xi(3)*(xi(1)+1)*(xi(3)+1))/2;
        dNdxi3 = -(xi(1)*(2*xi(3)+1)*(xi(2)^2-1)*(xi(1)+1))/4;
      case 19
        dNdxi1 = -(xi(1)*xi(2)*xi(3)*(xi(2)+1)*(xi(3)+1))/2;
        dNdxi2 = -(xi(3)*(2*xi(2)+1)*(xi(1)^2-1)*(xi(3)+1))/4;
        dNdxi3 = -(xi(2)*(2*xi(3)+1)*(xi(1)^2-1)*(xi(2)+1))/4;
      case 20
        dNdxi1 = -(xi(3)*(2*xi(1)-1)*(xi(2)^2-1)*(xi(3)+1))/4;
        dNdxi2 = -(xi(1)*xi(2)*xi(3)*(xi(1)-1)*(xi(3)+1))/2;
        dNdxi3 = -(xi(1)*(2*xi(3)+1)*(xi(2)^2-1)*(xi(1)-1))/4;
      case 21
        dNdxi1 = xi(1)*xi(3)*(xi(2)^2-1)*(xi(3)-1);
        dNdxi2 = xi(2)*xi(3)*(xi(1)^2-1)*(xi(3)-1);
        dNdxi3 = ((2*xi(3)-1)*(xi(1)^2-1)*(xi(2)^2-1))/2;
      case 22
        dNdxi1 = xi(1)*xi(2)*(xi(3)^2-1)*(xi(2)-1);
        dNdxi2 = ((2*xi(2)-1)*(xi(1)^2-1)*(xi(3)^2-1))/2;
        dNdxi3 = xi(2)*xi(3)*(xi(1)^2-1)*(xi(2)-1);
      case 23
        dNdxi1 = ((2*xi(1)+1)*(xi(2)^2-1)*(xi(3)^2-1))/2;
        dNdxi2 = xi(1)*xi(2)*(xi(3)^2-1)*(xi(1)+1);
        dNdxi3 = xi(1)*xi(3)*(xi(2)^2-1)*(xi(1)+1);
      case 24
        dNdxi1 = xi(1)*xi(2)*(xi(3)^2-1)*(xi(2)+1);
        dNdxi2 = ((2*xi(2)+1)*(xi(1)^2-1)*(xi(3)^2-1))/2;
        dNdxi3 = xi(2)*xi(3)*(xi(1)^2-1)*(xi(2)+1);
      case 25
        dNdxi1 = ((2*xi(1)-1)*(xi(2)^2-1)*(xi(3)^2-1))/2;
        dNdxi2 = xi(1)*xi(2)*(xi(3)^2-1)*(xi(1)-1);
        dNdxi3 = xi(1)*xi(3)*(xi(2)^2-1)*(xi(1)-1);
      case 26
        dNdxi1 = xi(1)*xi(3)*(xi(2)^2-1)*(xi(3)+1);
        dNdxi2 = xi(2)*xi(3)*(xi(1)^2-1)*(xi(3)+1);
        dNdxi3 = ((2*xi(3)+1)*(xi(1)^2-1)*(xi(2)^2-1))/2;
      case 27
        dNdxi1 = -2*xi(1)*(xi(2)^2-1)*(xi(3)^2-1);
        dNdxi2 = -2*xi(2)*(xi(1)^2-1)*(xi(3)^2-1);
        dNdxi3 = -2*xi(3)*(xi(1)^2-1)*(xi(2)^2-1);
    end

    if     ( i_eval==2 )   % x-derivative.

      vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:,2)*dNdxi2 + aInvJac(:,3)*dNdxi3;

    elseif ( i_eval==3 )   % y-derivative.

      vBase = aInvJac(:,4)*dNdxi1 + aInvJac(:,5)*dNdxi2 + aInvJac(:,6)*dNdxi3;

    elseif ( i_eval==4 )   % z-derivative.

      vBase = aInvJac(:,7)*dNdxi1 + aInvJac(:,8)*dNdxi2 + aInvJac(:,9)*dNdxi3;
    end

  case {22,23,24,32,33,34,42,43,44}   % Evaluation of second order derivatives.

    if( any(any(abs([aInvJac(:,[2 3 4 6 7 8])])>eps*1e2)) )
      warning('sf_hex_Q2: 2nd derivatives for non-rectangular cells shapes not supported.')
    end

    switch i_dof

      case 1
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) - 1)*(xi(3) - 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3) - 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2) - 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3) - 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) - 1)*(xi(3) - 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1) - 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2) - 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1) - 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) - 1)*(xi(2) - 1))/4;

      case 2
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) - 1)*(xi(3) - 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3) - 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2) - 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3) - 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) + 1)*(xi(3) - 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1) + 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2) - 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1) + 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) + 1)*(xi(2) - 1))/4;

      case 3
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) + 1)*(xi(3) - 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3) - 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2) + 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3) - 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) + 1)*(xi(3) - 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1) + 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2) + 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1) + 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) + 1)*(xi(2) + 1))/4;

      case 4
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) + 1)*(xi(3) - 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3) - 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2) + 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3) - 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) - 1)*(xi(3) - 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1) - 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2) + 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1) - 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) - 1)*(xi(2) + 1))/4;

      case 5
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) - 1)*(xi(3) + 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3) + 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2) - 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3) + 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) - 1)*(xi(3) + 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1) - 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2) - 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1) - 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) - 1)*(xi(2) - 1))/4;

      case 6
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) - 1)*(xi(3) + 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3) + 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2) - 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3) + 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) + 1)*(xi(3) + 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1) + 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2) - 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1) + 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) + 1)*(xi(2) - 1))/4;

      case 7
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) + 1)*(xi(3) + 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3) + 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2) + 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3) + 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) + 1)*(xi(3) + 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1) + 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2) + 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1) + 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) + 1)*(xi(2) + 1))/4;

      case 8
        d2Ndxi1dxi1 = (xi(2)*xi(3)*(xi(2) + 1)*(xi(3) + 1))/4;
        d2Ndxi2dxi1 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3) + 1))/8;
        d2Ndxi3dxi1 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2) + 1))/8;
        d2Ndxi1dxi2 = (xi(3)*(2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3) + 1))/8;
        d2Ndxi2dxi2 = (xi(1)*xi(3)*(xi(1) - 1)*(xi(3) + 1))/4;
        d2Ndxi3dxi2 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1) - 1))/8;
        d2Ndxi1dxi3 = (xi(2)*(2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2) + 1))/8;
        d2Ndxi2dxi3 = (xi(1)*(2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1) - 1))/8;
        d2Ndxi3dxi3 = (xi(1)*xi(2)*(xi(1) - 1)*(xi(2) + 1))/4;

      case 9
        d2Ndxi1dxi1 = -(xi(2)*xi(3)*(xi(2) - 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi1 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi1 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(2) - 1))/2;
        d2Ndxi1dxi2 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi2 = -(xi(3)*(xi(1)^2 - 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi2 = -((2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1)^2 - 1))/4;
        d2Ndxi1dxi3 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi3 = -((2*xi(2) - 1)*(2*xi(3) - 1)*(xi(1)^2 - 1))/4;
        d2Ndxi3dxi3 = -(xi(2)*(xi(1)^2 - 1)*(xi(2) - 1))/2;

      case 10
        d2Ndxi1dxi1 = -(xi(3)*(xi(2)^2 - 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi1 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi1 = -((2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2)^2 - 1))/4;
        d2Ndxi1dxi2 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi2 = -(xi(1)*xi(3)*(xi(1) + 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(1) + 1))/2;
        d2Ndxi1dxi3 = -((2*xi(1) + 1)*(2*xi(3) - 1)*(xi(2)^2 - 1))/4;
        d2Ndxi2dxi3 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*(xi(2)^2 - 1)*(xi(1) + 1))/2;

      case 11
        d2Ndxi1dxi1 = -(xi(2)*xi(3)*(xi(2) + 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi1 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi1 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(2) + 1))/2;
        d2Ndxi1dxi2 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi2 = -(xi(3)*(xi(1)^2 - 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi2 = -((2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1)^2 - 1))/4;
        d2Ndxi1dxi3 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi3 = -((2*xi(2) + 1)*(2*xi(3) - 1)*(xi(1)^2 - 1))/4;
        d2Ndxi3dxi3 = -(xi(2)*(xi(1)^2 - 1)*(xi(2) + 1))/2;

      case 12
        d2Ndxi1dxi1 = -(xi(3)*(xi(2)^2 - 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi1 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi1 = -((2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2)^2 - 1))/4;
        d2Ndxi1dxi2 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(3) - 1))/2;
        d2Ndxi2dxi2 = -(xi(1)*xi(3)*(xi(1) - 1)*(xi(3) - 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(1) - 1))/2;
        d2Ndxi1dxi3 = -((2*xi(1) - 1)*(2*xi(3) - 1)*(xi(2)^2 - 1))/4;
        d2Ndxi2dxi3 = -(xi(1)*xi(2)*(2*xi(3) - 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*(xi(2)^2 - 1)*(xi(1) - 1))/2;

      case 13
        d2Ndxi1dxi1 = -(xi(2)*(xi(3)^2 - 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi1 = -((2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3)^2 - 1))/4;
        d2Ndxi3dxi1 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(2) - 1))/2;
        d2Ndxi1dxi2 = -((2*xi(1) - 1)*(2*xi(2) - 1)*(xi(3)^2 - 1))/4;
        d2Ndxi2dxi2 = -(xi(1)*(xi(3)^2 - 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(1) - 1))/2;
        d2Ndxi1dxi3 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi3 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*xi(2)*(xi(1) - 1)*(xi(2) - 1))/2;

      case 14
        d2Ndxi1dxi1 = -(xi(2)*(xi(3)^2 - 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi1 = -((2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3)^2 - 1))/4;
        d2Ndxi3dxi1 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(2) - 1))/2;
        d2Ndxi1dxi2 = -((2*xi(1) + 1)*(2*xi(2) - 1)*(xi(3)^2 - 1))/4;
        d2Ndxi2dxi2 = -(xi(1)*(xi(3)^2 - 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(1) + 1))/2;
        d2Ndxi1dxi3 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi3 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*xi(2)*(xi(1) + 1)*(xi(2) - 1))/2;

      case 15
        d2Ndxi1dxi1 = -(xi(2)*(xi(3)^2 - 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi1 = -((2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3)^2 - 1))/4;
        d2Ndxi3dxi1 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(2) + 1))/2;
        d2Ndxi1dxi2 = -((2*xi(1) + 1)*(2*xi(2) + 1)*(xi(3)^2 - 1))/4;
        d2Ndxi2dxi2 = -(xi(1)*(xi(3)^2 - 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(1) + 1))/2;
        d2Ndxi1dxi3 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi3 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*xi(2)*(xi(1) + 1)*(xi(2) + 1))/2;

      case 16
        d2Ndxi1dxi1 = -(xi(2)*(xi(3)^2 - 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi1 = -((2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3)^2 - 1))/4;
        d2Ndxi3dxi1 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(2) + 1))/2;
        d2Ndxi1dxi2 = -((2*xi(1) - 1)*(2*xi(2) + 1)*(xi(3)^2 - 1))/4;
        d2Ndxi2dxi2 = -(xi(1)*(xi(3)^2 - 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(1) - 1))/2;
        d2Ndxi1dxi3 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi3 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*xi(2)*(xi(1) - 1)*(xi(2) + 1))/2;

      case 17
        d2Ndxi1dxi1 = -(xi(2)*xi(3)*(xi(2) - 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi1 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi1 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(2) - 1))/2;
        d2Ndxi1dxi2 = -(xi(1)*xi(3)*(2*xi(2) - 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi2 = -(xi(3)*(xi(1)^2 - 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi2 = -((2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1)^2 - 1))/4;
        d2Ndxi1dxi3 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(2) - 1))/2;
        d2Ndxi2dxi3 = -((2*xi(2) - 1)*(2*xi(3) + 1)*(xi(1)^2 - 1))/4;
        d2Ndxi3dxi3 = -(xi(2)*(xi(1)^2 - 1)*(xi(2) - 1))/2;

      case 18
        d2Ndxi1dxi1 = -(xi(3)*(xi(2)^2 - 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi1 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi1 = -((2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2)^2 - 1))/4;
        d2Ndxi1dxi2 = -(xi(2)*xi(3)*(2*xi(1) + 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi2 = -(xi(1)*xi(3)*(xi(1) + 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(1) + 1))/2;
        d2Ndxi1dxi3 = -((2*xi(1) + 1)*(2*xi(3) + 1)*(xi(2)^2 - 1))/4;
        d2Ndxi2dxi3 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(1) + 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*(xi(2)^2 - 1)*(xi(1) + 1))/2;

      case 19
        d2Ndxi1dxi1 = -(xi(2)*xi(3)*(xi(2) + 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi1 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi1 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(2) + 1))/2;
        d2Ndxi1dxi2 = -(xi(1)*xi(3)*(2*xi(2) + 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi2 = -(xi(3)*(xi(1)^2 - 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi2 = -((2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1)^2 - 1))/4;
        d2Ndxi1dxi3 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(2) + 1))/2;
        d2Ndxi2dxi3 = -((2*xi(2) + 1)*(2*xi(3) + 1)*(xi(1)^2 - 1))/4;
        d2Ndxi3dxi3 = -(xi(2)*(xi(1)^2 - 1)*(xi(2) + 1))/2;

      case 20
        d2Ndxi1dxi1 = -(xi(3)*(xi(2)^2 - 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi1 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi1 = -((2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2)^2 - 1))/4;
        d2Ndxi1dxi2 = -(xi(2)*xi(3)*(2*xi(1) - 1)*(xi(3) + 1))/2;
        d2Ndxi2dxi2 = -(xi(1)*xi(3)*(xi(1) - 1)*(xi(3) + 1))/2;
        d2Ndxi3dxi2 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(1) - 1))/2;
        d2Ndxi1dxi3 = -((2*xi(1) - 1)*(2*xi(3) + 1)*(xi(2)^2 - 1))/4;
        d2Ndxi2dxi3 = -(xi(1)*xi(2)*(2*xi(3) + 1)*(xi(1) - 1))/2;
        d2Ndxi3dxi3 = -(xi(1)*(xi(2)^2 - 1)*(xi(1) - 1))/2;

      case 21
        d2Ndxi1dxi1 = xi(3)*(xi(2)^2 - 1)*(xi(3) - 1);
        d2Ndxi2dxi1 = 2*xi(1)*xi(2)*xi(3)*(xi(3) - 1);
        d2Ndxi3dxi1 = xi(1)*(2*xi(3) - 1)*(xi(2)^2 - 1);
        d2Ndxi1dxi2 = 2*xi(1)*xi(2)*xi(3)*(xi(3) - 1);
        d2Ndxi2dxi2 = xi(3)*(xi(1)^2 - 1)*(xi(3) - 1);
        d2Ndxi3dxi2 = xi(2)*(2*xi(3) - 1)*(xi(1)^2 - 1);
        d2Ndxi1dxi3 = xi(1)*(2*xi(3) - 1)*(xi(2)^2 - 1);
        d2Ndxi2dxi3 = xi(2)*(2*xi(3) - 1)*(xi(1)^2 - 1);
        d2Ndxi3dxi3 = (xi(1)^2 - 1)*(xi(2)^2 - 1);

      case 22
        d2Ndxi1dxi1 = xi(2)*(xi(3)^2 - 1)*(xi(2) - 1);
        d2Ndxi2dxi1 = xi(1)*(2*xi(2) - 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi1 = 2*xi(1)*xi(2)*xi(3)*(xi(2) - 1);
        d2Ndxi1dxi2 = xi(1)*(2*xi(2) - 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi2 = (xi(1)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi2 = xi(3)*(2*xi(2) - 1)*(xi(1)^2 - 1);
        d2Ndxi1dxi3 = 2*xi(1)*xi(2)*xi(3)*(xi(2) - 1);
        d2Ndxi2dxi3 = xi(3)*(2*xi(2) - 1)*(xi(1)^2 - 1);
        d2Ndxi3dxi3 = xi(2)*(xi(1)^2 - 1)*(xi(2) - 1);

      case 23
        d2Ndxi1dxi1 = (xi(2)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi1 = xi(2)*(2*xi(1) + 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi1 = xi(3)*(2*xi(1) + 1)*(xi(2)^2 - 1);
        d2Ndxi1dxi2 = xi(2)*(2*xi(1) + 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi2 = xi(1)*(xi(3)^2 - 1)*(xi(1) + 1);
        d2Ndxi3dxi2 = 2*xi(1)*xi(2)*xi(3)*(xi(1) + 1);
        d2Ndxi1dxi3 = xi(3)*(2*xi(1) + 1)*(xi(2)^2 - 1);
        d2Ndxi2dxi3 = 2*xi(1)*xi(2)*xi(3)*(xi(1) + 1);
        d2Ndxi3dxi3 = xi(1)*(xi(2)^2 - 1)*(xi(1) + 1);

      case 24
        d2Ndxi1dxi1 = xi(2)*(xi(3)^2 - 1)*(xi(2) + 1);
        d2Ndxi2dxi1 = xi(1)*(2*xi(2) + 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi1 = 2*xi(1)*xi(2)*xi(3)*(xi(2) + 1);
        d2Ndxi1dxi2 = xi(1)*(2*xi(2) + 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi2 = (xi(1)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi2 = xi(3)*(2*xi(2) + 1)*(xi(1)^2 - 1);
        d2Ndxi1dxi3 = 2*xi(1)*xi(2)*xi(3)*(xi(2) + 1);
        d2Ndxi2dxi3 = xi(3)*(2*xi(2) + 1)*(xi(1)^2 - 1);
        d2Ndxi3dxi3 = xi(2)*(xi(1)^2 - 1)*(xi(2) + 1);

      case 25
        d2Ndxi1dxi1 = (xi(2)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi1 = xi(2)*(2*xi(1) - 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi1 = xi(3)*(2*xi(1) - 1)*(xi(2)^2 - 1);
        d2Ndxi1dxi2 = xi(2)*(2*xi(1) - 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi2 = xi(1)*(xi(3)^2 - 1)*(xi(1) - 1);
        d2Ndxi3dxi2 = 2*xi(1)*xi(2)*xi(3)*(xi(1) - 1);
        d2Ndxi1dxi3 = xi(3)*(2*xi(1) - 1)*(xi(2)^2 - 1);
        d2Ndxi2dxi3 = 2*xi(1)*xi(2)*xi(3)*(xi(1) - 1);
        d2Ndxi3dxi3 = xi(1)*(xi(2)^2 - 1)*(xi(1) - 1);

      case 26
        d2Ndxi1dxi1 = xi(3)*(xi(2)^2 - 1)*(xi(3) + 1);
        d2Ndxi2dxi1 = 2*xi(1)*xi(2)*xi(3)*(xi(3) + 1);
        d2Ndxi3dxi1 = xi(1)*(2*xi(3) + 1)*(xi(2)^2 - 1);
        d2Ndxi1dxi2 = 2*xi(1)*xi(2)*xi(3)*(xi(3) + 1);
        d2Ndxi2dxi2 = xi(3)*(xi(1)^2 - 1)*(xi(3) + 1);
        d2Ndxi3dxi2 = xi(2)*(2*xi(3) + 1)*(xi(1)^2 - 1);
        d2Ndxi1dxi3 = xi(1)*(2*xi(3) + 1)*(xi(2)^2 - 1);
        d2Ndxi2dxi3 = xi(2)*(2*xi(3) + 1)*(xi(1)^2 - 1);
        d2Ndxi3dxi3 = (xi(1)^2 - 1)*(xi(2)^2 - 1);

      case 27
        d2Ndxi1dxi1 = -2*(xi(2)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi2dxi1 = -4*xi(1)*xi(2)*(xi(3)^2 - 1);
        d2Ndxi3dxi1 = -4*xi(1)*xi(3)*(xi(2)^2 - 1);
        d2Ndxi1dxi2 = -4*xi(1)*xi(2)*(xi(3)^2 - 1);
        d2Ndxi2dxi2 = -2*(xi(1)^2 - 1)*(xi(3)^2 - 1);
        d2Ndxi3dxi2 = -4*xi(2)*xi(3)*(xi(1)^2 - 1);
        d2Ndxi1dxi3 = -4*xi(1)*xi(3)*(xi(2)^2 - 1);
        d2Ndxi2dxi3 = -4*xi(2)*xi(3)*(xi(1)^2 - 1);
        d2Ndxi3dxi3 = -2*(xi(1)^2 - 1)*(xi(2)^2 - 1);

    end

    switch( i_eval )
      case 22
        vBase = aInvJac(:,1).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:,2)*d2Ndxi2dxi1 + aInvJac(:,3)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,2).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:,2)*d2Ndxi2dxi2 + aInvJac(:,3)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,3).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:,2)*d2Ndxi2dxi3 + aInvJac(:,3)*d2Ndxi3dxi3 );
      case 33
        vBase = aInvJac(:,4).*( aInvJac(:,4)*d2Ndxi1dxi1 + aInvJac(:,5)*d2Ndxi2dxi1 + aInvJac(:,6)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,5).*( aInvJac(:,4)*d2Ndxi1dxi2 + aInvJac(:,5)*d2Ndxi2dxi2 + aInvJac(:,6)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,6).*( aInvJac(:,4)*d2Ndxi1dxi3 + aInvJac(:,5)*d2Ndxi2dxi3 + aInvJac(:,6)*d2Ndxi3dxi3 );

      case 44
        vBase = aInvJac(:,7).*( aInvJac(:,7)*d2Ndxi1dxi1 + aInvJac(:,8)*d2Ndxi2dxi1 + aInvJac(:,9)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,8).*( aInvJac(:,7)*d2Ndxi1dxi2 + aInvJac(:,8)*d2Ndxi2dxi2 + aInvJac(:,9)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,9).*( aInvJac(:,7)*d2Ndxi1dxi3 + aInvJac(:,8)*d2Ndxi2dxi3 + aInvJac(:,9)*d2Ndxi3dxi3 );

      case {23,32}
        vBase = aInvJac(:,4).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:,2)*d2Ndxi2dxi1 + aInvJac(:,3)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,5).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:,2)*d2Ndxi2dxi2 + aInvJac(:,3)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,6).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:,2)*d2Ndxi2dxi3 + aInvJac(:,3)*d2Ndxi3dxi3 );

      case {24,42}
        vBase = aInvJac(:,7).*( aInvJac(:,1)*d2Ndxi1dxi1 + aInvJac(:,2)*d2Ndxi2dxi1 + aInvJac(:,3)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,8).*( aInvJac(:,1)*d2Ndxi1dxi2 + aInvJac(:,2)*d2Ndxi2dxi2 + aInvJac(:,3)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,9).*( aInvJac(:,1)*d2Ndxi1dxi3 + aInvJac(:,2)*d2Ndxi2dxi3 + aInvJac(:,3)*d2Ndxi3dxi3 );

      case {34,43}
        vBase = aInvJac(:,7).*( aInvJac(:,4)*d2Ndxi1dxi1 + aInvJac(:,5)*d2Ndxi2dxi1 + aInvJac(:,6)*d2Ndxi3dxi1 ) + ...
                aInvJac(:,8).*( aInvJac(:,4)*d2Ndxi1dxi2 + aInvJac(:,5)*d2Ndxi2dxi2 + aInvJac(:,6)*d2Ndxi3dxi2 ) + ...
                aInvJac(:,9).*( aInvJac(:,4)*d2Ndxi1dxi3 + aInvJac(:,5)*d2Ndxi2dxi3 + aInvJac(:,6)*d2Ndxi3dxi3 );
    end

  otherwise
    vBase = 0;

end
