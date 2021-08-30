function [ vBase, nLDof, xLDof, sfun ] = sf_quad_Q2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_QUAD_Q2 Biquadratic conforming shape function for quadrilaterals (Q2).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_QUAD_Q2( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming biquadratic Q2 shape functions on quadrilaterals
%   with values defined in the nodes, edges, and cell center. XI is
%   [-1..1]^2 reference coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  2             Number of space dimensions
%       n_vert      scalar:  4             Number of vertices per cell
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
%   See also SFLAG2

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [4 4 0 1];
xLDof = [-1  1 1 -1  0 1 0 -1 0; ...
         -1 -1 1  1 -1 0 1  0 0];
sfun  = 'sf_quad_Q2';


switch i_eval   % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.
      case 1
        vBase =  (1-xi(1))*(1-xi(2))*xi(1)*xi(2)/4;
      case 2
        vBase = -(1+xi(1))*(1-xi(2))*xi(1)*xi(2)/4;
      case 3
        vBase =  (1+xi(1))*(1+xi(2))*xi(1)*xi(2)/4;
      case 4
        vBase = -(1-xi(1))*(1+xi(2))*xi(1)*xi(2)/4;
      case 5
        vBase = -(1-xi(1)*xi(1))*(1-xi(2))*xi(2)/2;
      case 6
        vBase =  (1+xi(1))*(1-xi(2)*xi(2))*xi(1)/2;
      case 7
        vBase =  (1-xi(1)*xi(1))*(1+xi(2))*xi(2)/2;
      case 8
        vBase = -(1-xi(1))*(1-xi(2)*xi(2))*xi(1)/2;
      case 9
        vBase =  (1-xi(1)*xi(1))*(1-xi(2)*xi(2));
    end

  case {2,3}   % Evaluation of first order derivatives.

    switch i_dof   % Basis function to evaluate.
      case 1
        dNdxi1 =  (1-2*xi(1))*(1-xi(2))*xi(2)/4;
        dNdxi2 =  (1-xi(1))*(1-2*xi(2))*xi(1)/4;
      case 2
        dNdxi1 = -(1+2*xi(1))*(1-xi(2))*xi(2)/4;
        dNdxi2 = -(1+xi(1))*(1-2*xi(2))*xi(1)/4;
      case 3
        dNdxi1 =  (1+2*xi(1))*(1+xi(2))*xi(2)/4;
        dNdxi2 =  (1+xi(1))*(1+2*xi(2))*xi(1)/4;
      case 4
        dNdxi1 = -(1-2*xi(1))*(1+xi(2))*xi(2)/4;
        dNdxi2 = -(1-xi(1))*(1+2*xi(2))*xi(1)/4;
      case 5
        dNdxi1 =  (1-xi(2))*xi(1)*xi(2);
        dNdxi2 = -(1-xi(1)*xi(1))*(1-2*xi(2))/2;
      case 6
        dNdxi1 =  (1+2*xi(1))*(1-xi(2)*xi(2))/2;
        dNdxi2 = -(1+xi(1))*xi(1)*xi(2);
      case 7
        dNdxi1 = -(1+xi(2))*xi(1)*xi(2);
        dNdxi2 =  (1-xi(1)*xi(1))*(1+2*xi(2))/2;
      case 8
        dNdxi1 = -(1-2*xi(1))*(1-xi(2)*xi(2))/2;
        dNdxi2 =  (1-xi(1))*xi(1)*xi(2);
      case 9
        dNdxi1 = -2*(1-xi(2)*xi(2))*xi(1);
        dNdxi2 = -2*(1-xi(1)*xi(1))*xi(2);
    end

    if     ( i_eval==2 )   % x-derivative.

      vBase = aInvJac(:,1)*dNdxi1+aInvJac(:,2)*dNdxi2;

    elseif ( i_eval==3 )   % y-derivative.

      vBase = aInvJac(:,3)*dNdxi1+aInvJac(:,4)*dNdxi2;
    end

  case {22,23,32,33}

    if( any(abs([aInvJac(:,2);aInvJac(:,3)])>eps*1e2) )
      warning('sf_quad_Q2: 2nd derivatives for non-rectangular cells shapes not supported.')
    end

    switch i_dof   % Basis function to evaluate.
      case 1
        d2Ndxi12    = xi(2)^2/2 - xi(2)/2;
        d2Ndxi1dxi2 = ((2*xi(1) - 1)*(2*xi(2) - 1))/4;
        d2Ndxi2dxi1 = ((2*xi(1) - 1)*(2*xi(2) - 1))/4;
        d2Ndxi22    = xi(1)^2/2 - xi(1)/2;
      case 2
        d2Ndxi12    = xi(2)^2/2 - xi(2)/2;
        d2Ndxi1dxi2 = ((2*xi(1) + 1)*(2*xi(2) - 1))/4;
        d2Ndxi2dxi1 = ((2*xi(1) + 1)*(2*xi(2) - 1))/4;
        d2Ndxi22    = xi(1)^2/2 + xi(1)/2;
      case 3
        d2Ndxi12    = xi(2)^2/2 + xi(2)/2;
        d2Ndxi1dxi2 = ((2*xi(1) + 1)*(2*xi(2) + 1))/4;
        d2Ndxi2dxi1 = ((2*xi(1) + 1)*(2*xi(2) + 1))/4;
        d2Ndxi22    = xi(1)^2/2 + xi(1)/2;
      case 4
        d2Ndxi12    = xi(2)^2/2 + xi(2)/2;
        d2Ndxi1dxi2 = ((2*xi(1) - 1)*(2*xi(2) + 1))/4;
        d2Ndxi2dxi1 = ((2*xi(1) - 1)*(2*xi(2) + 1))/4;
        d2Ndxi22    = xi(1)^2/2 - xi(1)/2;
      case 5
        d2Ndxi12    = xi(2) - xi(2)^2;
        d2Ndxi1dxi2 = xi(1) - 2*xi(1)*xi(2);
        d2Ndxi2dxi1 = xi(1) - 2*xi(1)*xi(2);
        d2Ndxi22    = 1 - xi(1)^2;
      case 6
        d2Ndxi12    = 1 - xi(2)^2;
        d2Ndxi1dxi2 = -xi(2)*(2*xi(1) + 1);
        d2Ndxi2dxi1 = -xi(2)*(2*xi(1) + 1);
        d2Ndxi22    = - xi(1)^2 - xi(1);
      case 7
        d2Ndxi12    = - xi(2)^2 - xi(2);
        d2Ndxi1dxi2 = -xi(1)*(2*xi(2) + 1);
        d2Ndxi2dxi1 = -xi(1)*(2*xi(2) + 1);
        d2Ndxi22    = 1 - xi(1)^2;
      case 8
        d2Ndxi12    = 1 - xi(2)^2;
        d2Ndxi1dxi2 = xi(2) - 2*xi(1)*xi(2);
        d2Ndxi2dxi1 = xi(2) - 2*xi(1)*xi(2);
        d2Ndxi22    = xi(1) - xi(1)^2;
      case 9
        d2Ndxi12    = 2*xi(2)^2 - 2;
        d2Ndxi1dxi2 = 4*xi(1)*xi(2);
        d2Ndxi2dxi1 = 4*xi(1)*xi(2);
        d2Ndxi22    = 2*xi(1)^2 - 2;
    end

    if ( i_eval==22 )   % xx-derivative.

      vBase  = aInvJac(:,1).*( aInvJac(:,1)*d2Ndxi12    + aInvJac(:,2)*d2Ndxi1dxi2 ) + ...
               aInvJac(:,2).*( aInvJac(:,1)*d2Ndxi2dxi1 + aInvJac(:,2)*d2Ndxi22    );

    elseif ( i_eval==23 )   % xy-derivative.

      vBase  = aInvJac(:,3).*( aInvJac(:,1)*d2Ndxi12    + aInvJac(:,2)*d2Ndxi1dxi2 ) + ...
               aInvJac(:,4).*( aInvJac(:,1)*d2Ndxi2dxi1 + aInvJac(:,2)*d2Ndxi22    );

    elseif ( i_eval==32 )   % yx-derivative.

      vBase  = aInvJac(:,1).*( aInvJac(:,3)*d2Ndxi12    + aInvJac(:,4)*d2Ndxi1dxi2 ) + ...
               aInvJac(:,2).*( aInvJac(:,3)*d2Ndxi2dxi1 + aInvJac(:,4)*d2Ndxi22    );

    elseif ( i_eval==33 )   % yy-derivative.

      vBase  = aInvJac(:,3).*( aInvJac(:,3)*d2Ndxi12    + aInvJac(:,4)*d2Ndxi1dxi2 ) + ...
               aInvJac(:,4).*( aInvJac(:,3)*d2Ndxi2dxi1 + aInvJac(:,4)*d2Ndxi22    );
    end

  otherwise
    vBase = 0;

end
