function [ vBase, nLDof, xLDof, sfun ] = sf_hex_Q1( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_HEX_Q1 Trilinear conforming shape function for hexahedrons (Q1).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_HEX_Q1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming trilinear Q1 shape functions on hexahedrons with
%   values defined in the nodes. XI is [-1..1]^3 reference coordinates.
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
%   See also SFLAG1, SF_HEX_Q1NC

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [8 0 0 0];
xLDof = [-1  1  1 -1 -1  1  1 -1; ...
         -1 -1  1  1 -1 -1  1  1; ...
         -1 -1 -1 -1  1  1  1  1];
sfun  = 'sf_hex_Q1';


switch i_eval   % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = (1-xi(1))*(1-xi(2))*(1-xi(3))/8;
      case 2
        vBase = (1+xi(1))*(1-xi(2))*(1-xi(3))/8;
      case 3
        vBase = (1+xi(1))*(1+xi(2))*(1-xi(3))/8;
      case 4
        vBase = (1-xi(1))*(1+xi(2))*(1-xi(3))/8;
      case 5
        vBase = (1-xi(1))*(1-xi(2))*(1+xi(3))/8;
      case 6
        vBase = (1+xi(1))*(1-xi(2))*(1+xi(3))/8;
      case 7
        vBase = (1+xi(1))*(1+xi(2))*(1+xi(3))/8;
      case 8
        vBase = (1-xi(1))*(1+xi(2))*(1+xi(3))/8;
    end

  case {2,3,4}   % Evaluation of first order derivatives.

    switch i_dof   % Basis function to evaluate.

      case 1
        dNdxi1 = -(1-xi(2))*(1-xi(3))/8;
        dNdxi2 = -(1-xi(1))*(1-xi(3))/8;
        dNdxi3 = -(1-xi(1))*(1-xi(2))/8;
      case 2
        dNdxi1 =  (1-xi(2))*(1-xi(3))/8;
        dNdxi2 = -(1+xi(1))*(1-xi(3))/8;
        dNdxi3 = -(1+xi(1))*(1-xi(2))/8;
      case 3
        dNdxi1 =  (1+xi(2))*(1-xi(3))/8;
        dNdxi2 =  (1+xi(1))*(1-xi(3))/8;
        dNdxi3 = -(1+xi(1))*(1+xi(2))/8;
      case 4
        dNdxi1 = -(1+xi(2))*(1-xi(3))/8;
        dNdxi2 =  (1-xi(1))*(1-xi(3))/8;
        dNdxi3 = -(1-xi(1))*(1+xi(2))/8;
      case 5
        dNdxi1 = -(1-xi(2))*(1+xi(3))/8;
        dNdxi2 = -(1-xi(1))*(1+xi(3))/8;
        dNdxi3 =  (1-xi(1))*(1-xi(2))/8;
      case 6
        dNdxi1 =  (1-xi(2))*(1+xi(3))/8;
        dNdxi2 = -(1+xi(1))*(1+xi(3))/8;
        dNdxi3 =  (1+xi(1))*(1-xi(2))/8;
      case 7
        dNdxi1 =  (1+xi(2))*(1+xi(3))/8;
        dNdxi2 =  (1+xi(1))*(1+xi(3))/8;
        dNdxi3 =  (1+xi(1))*(1+xi(2))/8;
      case 8
        dNdxi1 = -(1+xi(2))*(1+xi(3))/8;
        dNdxi2 =  (1-xi(1))*(1+xi(3))/8;
        dNdxi3 =  (1-xi(1))*(1+xi(2))/8;
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
      warning('sf_hex_Q1: 2nd derivatives for non-rectangular cells shapes not supported.')
    end

    switch i_dof

      case 1
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 1/8 - xi(3)/8;
        d2Ndxi3dxi1 = 1/8 - xi(2)/8;
        d2Ndxi1dxi2 = 1/8 - xi(3)/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 1/8 - xi(1)/8;
        d2Ndxi1dxi3 = 1/8 - xi(2)/8;
        d2Ndxi2dxi3 = 1/8 - xi(1)/8;
        d2Ndxi3dxi3 = 0;

      case 2
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = xi(3)/8 - 1/8;
        d2Ndxi3dxi1 = xi(2)/8 - 1/8;
        d2Ndxi1dxi2 = xi(3)/8 - 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = xi(1)/8 + 1/8;
        d2Ndxi1dxi3 = xi(2)/8 - 1/8;
        d2Ndxi2dxi3 = xi(1)/8 + 1/8;
        d2Ndxi3dxi3 = 0;

      case 3
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = 1/8 - xi(3)/8;
        d2Ndxi3dxi1 = - xi(2)/8 - 1/8;
        d2Ndxi1dxi2 = 1/8 - xi(3)/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = - xi(1)/8 - 1/8;
        d2Ndxi1dxi3 = - xi(2)/8 - 1/8;
        d2Ndxi2dxi3 = - xi(1)/8 - 1/8;
        d2Ndxi3dxi3 = 0;

      case 4
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = xi(3)/8 - 1/8;
        d2Ndxi3dxi1 = xi(2)/8 + 1/8;
        d2Ndxi1dxi2 = xi(3)/8 - 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = xi(1)/8 - 1/8;
        d2Ndxi1dxi3 = xi(2)/8 + 1/8;
        d2Ndxi2dxi3 = xi(1)/8 - 1/8;
        d2Ndxi3dxi3 = 0;

      case 5
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = xi(3)/8 + 1/8;
        d2Ndxi3dxi1 = xi(2)/8 - 1/8;
        d2Ndxi1dxi2 = xi(3)/8 + 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = xi(1)/8 - 1/8;
        d2Ndxi1dxi3 = xi(2)/8 - 1/8;
        d2Ndxi2dxi3 = xi(1)/8 - 1/8;
        d2Ndxi3dxi3 = 0;

      case 6
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = - xi(3)/8 - 1/8;
        d2Ndxi3dxi1 = 1/8 - xi(2)/8;
        d2Ndxi1dxi2 = - xi(3)/8 - 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = - xi(1)/8 - 1/8;
        d2Ndxi1dxi3 = 1/8 - xi(2)/8;
        d2Ndxi2dxi3 = - xi(1)/8 - 1/8;
        d2Ndxi3dxi3 = 0;

      case 7
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = xi(3)/8 + 1/8;
        d2Ndxi3dxi1 = xi(2)/8 + 1/8;
        d2Ndxi1dxi2 = xi(3)/8 + 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = xi(1)/8 + 1/8;
        d2Ndxi1dxi3 = xi(2)/8 + 1/8;
        d2Ndxi2dxi3 = xi(1)/8 + 1/8;
        d2Ndxi3dxi3 = 0;

      case 8
        d2Ndxi1dxi1 = 0;
        d2Ndxi2dxi1 = - xi(3)/8 - 1/8;
        d2Ndxi3dxi1 = - xi(2)/8 - 1/8;
        d2Ndxi1dxi2 = - xi(3)/8 - 1/8;
        d2Ndxi2dxi2 = 0;
        d2Ndxi3dxi2 = 1/8 - xi(1)/8;
        d2Ndxi1dxi3 = - xi(2)/8 - 1/8;
        d2Ndxi2dxi3 = 1/8 - xi(1)/8;
        d2Ndxi3dxi3 = 0;

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
