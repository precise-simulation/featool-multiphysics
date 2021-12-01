function [ vBase, nLDof, xLDof, sfun ] = sf_simp_P1bub( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_SIMP_P1BUB Linear Lagrange shape function for simplices with bubble (P1+).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_P1BUB( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming linear P1 Lagrange shape functions on simplices
%   an additional with bubble function. XI Barycentric coordinates.
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
%   See also SF_SIMP_P1

% Copyright 2013-2021 Precise Simulation, Ltd.


sfun = 'sf_simp_P1bub';
[~,nLDof,xLDof] = sf_simp_P1( 0, n_sdim, n_vert );

nLDof(4) = 1;
switch n_sdim
  case 1
    xLDof = [ xLDof [1;1]/2 ];
  case 2
    xLDof = [ xLDof [1;1;1]/3 ];
  case 3
    xLDof = [ xLDof [1;1;1;1]/4 ];
end


% Evaluation type flag.
if( i_eval==1 )    % Evaluation of function values.

  if( n_sdim==1 )
    vBase = sf_simp_P2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase );

  elseif( n_sdim==2 )

    switch( i_dof )
      case 1
        vBase = (9*xi(2)*xi(3) - 1)*(xi(2) + xi(3) - 1);

      case 2
        vBase = xi(2) + 9*xi(2)*xi(3)*(xi(2) + xi(3) - 1);

      case 3
        vBase = xi(3) + 9*xi(2)*xi(3)*(xi(2) + xi(3) - 1);

      case 4
        vBase = -27*xi(2)*xi(3)*(xi(2) + xi(3) - 1);

    end

  else   % 3D.

    switch( i_dof )

      case 1
        vBase = (64*xi(2)*xi(3)*xi(4) - 1)*(xi(2) + xi(3) + xi(4) - 1);

      case 2
        vBase = xi(2)*(64*xi(3)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 1);

      case 3
        vBase = xi(3)*(64*xi(2)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 1);

      case 4
        vBase = xi(4)*(64*xi(2)*xi(3)*(xi(2) + xi(3) + xi(4) - 1) + 1);

      case 5
        vBase = -256*xi(2)*xi(3)*xi(4)*(xi(2) + xi(3) + xi(4) - 1);

    end

  end

elseif( i_eval>=2 && i_eval<=n_sdim+1 )   % Evaluation of first derivatives.

  if( n_sdim==1 )
    vBase = sf_simp_P2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase );

  elseif( n_sdim==2 )

    switch i_dof   % Basis function to evaluate.

      case 1
        dNdxi1 = 0;
        dNdxi2 = 18*xi(2)*xi(3) - 9*xi(3) + 9*xi(3)^2 - 1;
        dNdxi3 = 18*xi(2)*xi(3) - 9*xi(2) + 9*xi(2)^2 - 1;

      case 2
        dNdxi1 = 0;
        dNdxi2 = 18*xi(2)*xi(3) - 9*xi(3) + 9*xi(3)^2 + 1;
        dNdxi3 = 9*xi(2)*(xi(2) + 2*xi(3) - 1);

      case 3
        dNdxi1 = 0;
        dNdxi2 = 9*xi(3)*(2*xi(2) + xi(3) - 1);
        dNdxi3 = 18*xi(2)*xi(3) - 9*xi(2) + 9*xi(2)^2 + 1;

      case 4
        dNdxi1 = 0;
        dNdxi2 = -27*xi(3)*(2*xi(2) + xi(3) - 1);
        dNdxi3 = -27*xi(2)*(xi(2) + 2*xi(3) - 1);

    end

    if( i_eval==2 )

      vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:,2)*dNdxi2 + aInvJac(:,3)*dNdxi3;

    else

      vBase = aInvJac(:,4)*dNdxi1 + aInvJac(:,5)*dNdxi2 + aInvJac(:,6)*dNdxi3;

    end

  else   % 3D.

    switch i_dof   % Basis function to evaluate.

      case 1
        dNdxi1 = 0;
        dNdxi2 = 64*xi(3)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) - 1;
        dNdxi3 = 64*xi(2)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) - 1;
        dNdxi4 = 64*xi(2)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) - 1;

      case 2
        dNdxi1 = 0;
        dNdxi2 = 64*xi(3)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) + 1;
        dNdxi3 = 64*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);
        dNdxi4 = 64*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);

      case 3
        dNdxi1 = 0;
        dNdxi2 = 64*xi(3)*xi(4)*(2*xi(2) + xi(3) + xi(4) - 1);
        dNdxi3 = 64*xi(2)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) + 1;
        dNdxi4 = 64*xi(2)*xi(4)*(xi(2) + xi(3) + xi(4) - 1) + 64*xi(2)*xi(3)*xi(4) + 1;

      case 4
        dNdxi1 = 0;
        dNdxi2 = 64*xi(3)*xi(4)*(2*xi(2) + xi(3) + xi(4) - 1);
        dNdxi3 = 64*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);
        dNdxi4 = 64*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);

      case 5
        dNdxi1 = 0;
        dNdxi2 = -256*xi(3)*xi(4)*(2*xi(2) + xi(3) + xi(4) - 1);
        dNdxi3 = -256*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);
        dNdxi4 = -256*xi(2)*xi(4)*(xi(2) + 2*xi(3) + xi(4) - 1);

    end

    if( i_eval==2 )

      vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:,2)*dNdxi2 + aInvJac(:,3)*dNdxi3 + aInvJac(:,4)*dNdxi4;

    elseif( i_eval==3 )

      vBase = aInvJac(:,5)*dNdxi1 + aInvJac(:,6)*dNdxi2 + aInvJac(:,7)*dNdxi3 + aInvJac(:,8)*dNdxi4;

    else

      vBase = aInvJac(:,9)*dNdxi1 + aInvJac(:,10)*dNdxi2 + aInvJac(:,11)*dNdxi3 + aInvJac(:,12)*dNdxi4;

    end

  end

elseif( any(i_eval==[22 23 24 32 33 34 42 43 44]) )   % Evaluation of second derivatives.
  error('sf_simp_P1bub: second order derivative evaluation not supported.')

else

  vBase = 0;

end
