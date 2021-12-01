function [ vBase, nLDof, xLDof, sfun ] = sf_quad_Q1nc( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_QUAD_Q1NC Bilinear nonconforming shape function for quadrilaterals (Q1~).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_QUAD_Q1NC( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates nonconforming rotated bilinear Q1~ shape functions on quadrilaterals
%   with values defined on the edge midpoints. XI is [-1..1]^2 reference coordinates.
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
%   See also SFLAG1, SF_QUAD_Q1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [0 4 0 0];
xLDof = [ 0  1 0 -1; ...
         -1  0 1  0];
sfun  = 'sf_quad_Q1nc';


switch i_eval   % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = (-xi(1)^2+xi(2)^2-2*xi(2)+1)/4;
      case 2
        vBase = ( xi(1)^2-xi(2)^2+2*xi(1)+1)/4;
      case 3
        vBase = (-xi(1)^2+xi(2)^2+2*xi(2)+1)/4;
      case 4
        vBase = ( xi(1)^2-xi(2)^2-2*xi(1)+1)/4;
    end

  case {2,3}   % Evaluation of first order derivatives.

    switch i_dof   % Basis function to evaluate.

      case 1
        dNdxi1 = -xi(1)/2;
        dNdxi2 = (xi(2)-1)/2;
      case 2
        dNdxi1 = (xi(1)+1)/2;
        dNdxi2 = -xi(2)/2;
      case 3
        dNdxi1 = -xi(1)/2;
        dNdxi2 = (xi(2)+1)/2;
      case 4
        dNdxi1 = (xi(1)-1)/2;
        dNdxi2 = -xi(2)/2;
    end

    if     ( i_eval==2 )   % x-derivative.

      vBase = aInvJac(:,1)*dNdxi1+aInvJac(:,2)*dNdxi2;

    elseif ( i_eval==3 )   % y-derivative.

      vBase = aInvJac(:,3)*dNdxi1+aInvJac(:,4)*dNdxi2;
    end

  case {22,23,32,33}   % Evaluation of second order derivatives.

    if( any(abs([aInvJac(:,2);aInvJac(:,3)])>eps*1e2) )
      warning('sf_quad_Q1nc: 2nd derivatives for non-rectangular cells shapes not supported.')
    end

    switch i_dof   % Basis function to evaluate.

      case {1,3}
        d2Ndxi12    = -1/2;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi22    = 1/2;
      case {2,4}
        d2Ndxi12    = 1/2;
        d2Ndxi1dxi2 = 0;
        d2Ndxi2dxi1 = 0;
        d2Ndxi22    = -1/2;
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
