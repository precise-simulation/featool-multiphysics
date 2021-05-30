function [ vBase, nLDof, xLDof, sfun ] = sf_simp_P2bub( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_SIMP_P2BUB Quadratic Lagrange shape function for simplices with bubble (P2+).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_P2BUB( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming quadratic P2 Lagrange shape functions on simplices
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
%   See also SF_SIMP_P2

% Copyright 2013-2021 Precise Simulation, Ltd.


[~,nLDof,xLDof] = sf_simp_P2( 0, n_sdim, n_vert );

nLDof(4) = 1;
switch n_sdim
  case 2
    xLDof = [ xLDof [1;1;1]/3 ];
  case 3
    xLDof = [ xLDof [1;1;1;1]/4 ];
end
sfun = 'sf_simp_P2bub';


% Evaluation type flag.
vBase = 0;
if( nargin>=4 && i_dof==sum(nLDof) && n_sdim>1 )

  if( i_eval==1 )    % Evaluation of function values.

    switch n_sdim

      case 2

        vBase = 27*xi(1)*xi(2)*xi(3);

      case 3

        vBase = 256*xi(1)*xi(2)*xi(3)*xi(4);

    end

  elseif( i_eval>=2 && i_eval<=n_sdim+1 )   % Evaluation of first derivatives.

    switch n_sdim

      case 2

        dNdxi1 = 27*xi(2)*xi(3);
        dNdxi2 = 27*xi(1)*xi(3);
        dNdxi3 = 27*xi(1)*xi(2);

        if( i_eval==2 )
          vBase = aInvJac(:,1).*dNdxi1 + aInvJac(:,2).*dNdxi2 + aInvJac(:,3).*dNdxi3;
        else
          vBase = aInvJac(:,4).*dNdxi1 + aInvJac(:,5).*dNdxi2 + aInvJac(:,6).*dNdxi3;
        end

      case 3

        dNdxi1 = 256*xi(2)*xi(3)*xi(4);
        dNdxi2 = 256*xi(1)*xi(3)*xi(4);
        dNdxi3 = 256*xi(1)*xi(2)*xi(4);
        dNdxi4 = 256*xi(1)*xi(2)*xi(3);

        if( i_eval==2 )
          vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:, 2)*dNdxi2 + aInvJac(:, 3)*dNdxi3 + aInvJac(:, 4)*dNdxi4;
        elseif( i_eval==3 )
          vBase = aInvJac(:,5)*dNdxi1 + aInvJac(:, 6)*dNdxi2 + aInvJac(:, 7)*dNdxi3 + aInvJac(:, 8)*dNdxi4;
        else
          vBase = aInvJac(:,9)*dNdxi1 + aInvJac(:,10)*dNdxi2 + aInvJac(:,11)*dNdxi3 + aInvJac(:,12)*dNdxi4;
        end

    end

  elseif( any(i_eval==[22 23 24 32 33 34 42 43 44]) )   % Evaluation of second derivatives.

    switch n_sdim

      case 2

        d2Ndxi1dxi1 = 0;
        d2Ndxi1dxi2 = 27*xi(3);
        d2Ndxi1dxi3 = 27*xi(2);
        d2Ndxi2dxi1 = 27*xi(1);
        d2Ndxi2dxi2 = 0;
        d2Ndxi2dxi3 = 27*xi(3);
        d2Ndxi3dxi1 = 27*xi(1);
        d2Ndxi3dxi2 = 27*xi(2);
        d2Ndxi3dxi3 = 0;

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

      case 3

        d2Ndxi1dxi1 = 0;
        d2Ndxi1dxi2 = 256*xi(3)*xi(4);
        d2Ndxi1dxi3 = 256*xi(2)*xi(4);
        d2Ndxi1dxi4 = 256*xi(2)*xi(3);
        d2Ndxi2dxi1 = 256*xi(3)*xi(4);
        d2Ndxi2dxi2 = 0;
        d2Ndxi2dxi3 = 256*xi(1)*xi(4);
        d2Ndxi2dxi4 = 256*xi(1)*xi(3);
        d2Ndxi3dxi1 = 256*xi(2)*xi(4);
        d2Ndxi3dxi2 = 256*xi(1)*xi(4);
        d2Ndxi3dxi3 = 0;
        d2Ndxi3dxi4 = 256*xi(1)*xi(2);
        d2Ndxi4dxi1 = 256*xi(2)*xi(3);
        d2Ndxi4dxi2 = 256*xi(1)*xi(3);
        d2Ndxi4dxi3 = 256*xi(1)*xi(2);
        d2Ndxi4dxi4 = 0;

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

    end

  end

elseif( nargin>=4 )

  vBase = sf_simp_P2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase );

end
