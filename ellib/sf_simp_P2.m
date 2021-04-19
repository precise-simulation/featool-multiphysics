function [ vBase, nLDof, xLDof, sfun ] = sf_simp_P2( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_SIMP_P2 Quadratic Lagrange shape function for simplices (P2).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_P2( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming quadratic P2 Lagrange shape functions on simplices
%   with values defined in the nodes. XI Barycentric coordinates.
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
%   See also SFLAG2

% Copyright 2013-2021 Precise Simulation, Ltd.


switch n_sdim

  case 1
    nLDof = [ 2 0 0 1 ];
    xLDof = [ 1 0 0.5; ...
              0 1 0.5 ];
  case 2
    nLDof = [ 3 3 0 0 ];
    xLDof = [ 1 0 0 0.5 0   0.5; ...
              0 1 0 0.5 0.5 0  ; ...
              0 0 1 0   0.5 0.5 ];
  case 3
    nLDof = [ 4 6 0 0 ];
    xLDof = [ 1 0 0 0 0.5 0   0.5 0.5 0   0  ; ...
              0 1 0 0 0.5 0.5 0   0   0.5 0  ; ...
              0 0 1 0 0   0.5 0.5 0   0   0.5; ...
              0 0 0 1 0   0   0   0.5 0.5 0.5 ];
end
sfun = 'sf_simp_P2';


% Evaluation type flag.
if( i_eval==1 )    % Evaluation of function values.

    if( i_dof<=n_vert )   % Node dofs.

      vBase = xi(i_dof)*(2*xi(i_dof)-1);

    else   % Edge dofs.

      i = mod(i_dof-n_vert-1,min(n_vert,3))+1;
      j = mod(i,min(n_vert,3))+1;
      if ( i_dof>7 )
        j = 4;
      end

      vBase = 4*xi(i)*xi(j);

    end

elseif( i_eval>=2 && i_eval<=n_sdim+1 )   % Evaluation of first derivatives.

    if( i_dof<=n_vert )   % Node dofs.

      i_col = i_dof + (i_eval-2)*n_vert;
      dNdxi = 4*xi(i_dof)-1;

      vBase = aInvJac(:,i_col)*dNdxi;

    else   % Edge dofs.

      i = mod(i_dof-n_vert-1,min(n_vert,3))+1;
      j = mod(i,min(n_vert,3))+1;
      if ( i_dof>7 )
        j = 4;
      end
      i_col  = i + (i_eval-2)*n_vert;
      j_col  = j + (i_eval-2)*n_vert;
      dNdxii = 4*xi(j);
      dNdxij = 4*xi(i);

      vBase  = aInvJac(:,i_col)*dNdxii + aInvJac(:,j_col)*dNdxij;

    end

elseif( any(i_eval==[22 23 24 32 33 34 42 43 44]) )   % Evaluation of second derivatives.

    i_eval1 = floor(i_eval/10);
    i_eval2 = i_eval - 10*i_eval1;

    if( i_dof<=n_vert )   % Node dofs.

      d2Ndxi2 = 4;

      i_eval = floor(i_eval/10);
      i_col  = i_dof + (i_eval1-2)*n_vert;
      j_col  = i_dof + (i_eval2-2)*n_vert;
      vBase  = aInvJac(:,i_col).*aInvJac(:,j_col)*d2Ndxi2;

    else   % Edge dofs.

      d2Ndxii2    = 0;
      d2Ndxiidxij = 4;
      d2Ndxijdxii = d2Ndxiidxij;
      d2Ndxij2    = 0;

      i = mod(i_dof-n_vert-1,min(n_vert,3))+1;
      j = mod(i,min(n_vert,3))+1;
      if( i_dof>7 )
        j = 4;
      end

      i_col = i + (i_eval1-2)*n_vert;
      j_col = j + (i_eval1-2)*n_vert;

      k_col = i + (i_eval2-2)*n_vert;
      l_col = j + (i_eval2-2)*n_vert;

      vBase  = aInvJac(:,k_col).*( aInvJac(:,i_col)*d2Ndxii2    + aInvJac(:,j_col)*d2Ndxiidxij ) + ...
               aInvJac(:,l_col).*( aInvJac(:,i_col)*d2Ndxijdxii + aInvJac(:,j_col)*d2Ndxij2    );
    end

else

  vBase = 0;

end
