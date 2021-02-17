function [ vBase, nLDof, xLDof, sfun ] = sf_simp_N1( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_SIMP_N1 Linear vector (Nedelec) curl shape function for simplices.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_SIMP_N1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates linear vector Nedelec curl shape functions on simplices with
%   values defined in the nodes. XI is Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 2/3            Number of space dimensions
%       n_vert      scalar: 3/4            Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim+1]             Local coordinates of evaluation point
%       aInvJac     [n,n_sdim+1*n_sdim]    Inverse of transformation Jacobian
%       vBase       [n,1,2/3]              Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n,1,2/3]              Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [n_sdim,n_ldof]        Local coordinates of local dofs
%       sfun        string                 Function name of called shape function

% Copyright 2013-2021 Precise Simulation, Ltd.


sfun = 'sf_simp_N1';
if( n_sdim==1 )
  error( [sfun,': shape function not defined in 1D.'] )
end

if( n_sdim==2 )
  nLDof = [0 3 0 0];
  xLDof = [1/2 0   1/2;
           1/2 1/2 0  ;
           0   1/2 1/2];
else
  nLDof = [0 6 0 0];
  xLDof = [1/2 0   1/2 1/2 0   0  ;
           1/2 1/2 0   0   1/2 0  ;
           0   1/2 1/2 0   0   1/2;
           0   0   0   1/2 1/2 1/2];
end


switch i_eval    % Evaluation type flag.

  case 1         % Evaluation of function values.

    if( n_sdim==2 )
      j_dof = mod(i_dof,3) + 1;

      vBase = cat( 3, xi(i_dof)*aInvJac(:,j_dof)        - xi(j_dof)*aInvJac(:,i_dof), ...
                      xi(i_dof)*aInvJac(:,j_dof+n_vert) - xi(j_dof)*aInvJac(:,i_dof+n_vert) );

    else

      if( i_dof<=3 )
        j_dof = mod(i_dof,3) + 1;
      else
        i_dof = mod(i_dof-1,3) + 1;
        j_dof = 4;
      end

      vBase = cat( 3, xi(i_dof)*aInvJac(:,j_dof)          - xi(j_dof)*aInvJac(:,i_dof), ...
                      xi(i_dof)*aInvJac(:,j_dof+n_vert)   - xi(j_dof)*aInvJac(:,i_dof+n_vert), ...
                      xi(i_dof)*aInvJac(:,j_dof+2*n_vert) - xi(j_dof)*aInvJac(:,i_dof+2*n_vert) );

    end

  case {2,3,4}   % Evaluation of first derivatives.

    if( n_sdim==2 )

      j_dof = mod(i_dof,3) + 1;

      dNdxii1 =  aInvJac(:,j_dof);
      dNdxij1 = -aInvJac(:,i_dof);

      dNdxii2 =  aInvJac(:,j_dof+n_vert);
      dNdxij2 = -aInvJac(:,i_dof+n_vert);

      vBase = cat( 3, ( aInvJac(:,i_dof+3*(i_eval-2)).*dNdxii1 + aInvJac(:,j_dof+3*(i_eval-2)).*dNdxij1 ), ...
                      ( aInvJac(:,i_dof+3*(i_eval-2)).*dNdxii2 + aInvJac(:,j_dof+3*(i_eval-2)).*dNdxij2 ) );

    else

      if( i_dof<=3 )
        j_dof = mod(i_dof,3) + 1;
      else
        i_dof = mod(i_dof-1,3) + 1;
        j_dof = 4;
      end

      dNdxii1 =  aInvJac(:,j_dof);
      dNdxij1 = -aInvJac(:,i_dof);

      dNdxii2 =  aInvJac(:,j_dof+n_vert);
      dNdxij2 = -aInvJac(:,i_dof+n_vert);

      dNdxii3 =  aInvJac(:,j_dof+2*n_vert);
      dNdxij3 = -aInvJac(:,i_dof+2*n_vert);

      vBase = cat( 3, ( aInvJac(:,i_dof+4*(i_eval-2)).*dNdxii1 + aInvJac(:,j_dof+4*(i_eval-2)).*dNdxij1 ), ...
                      ( aInvJac(:,i_dof+4*(i_eval-2)).*dNdxii2 + aInvJac(:,j_dof+4*(i_eval-2)).*dNdxij2 ), ...;
                      ( aInvJac(:,i_dof+4*(i_eval-2)).*dNdxii3 + aInvJac(:,j_dof+4*(i_eval-2)).*dNdxij3 ) );

    end

  case {5}   % Evaluation of curl.

    if( n_sdim==2 )

      j_dof = mod(i_dof,3) + 1;

      dNdxii1 =  aInvJac(:,j_dof);
      dNdxij1 = -aInvJac(:,i_dof);

      dNdxii2 =  aInvJac(:,j_dof+n_vert);
      dNdxij2 = -aInvJac(:,i_dof+n_vert);

      vBase = ( aInvJac(:,i_dof)  .*dNdxii2 + aInvJac(:,j_dof)  .*dNdxij2 ) - ...
              ( aInvJac(:,i_dof+3).*dNdxii1 + aInvJac(:,j_dof+3).*dNdxij1 );

    else

      if( i_dof<=3 )
        j_dof = mod(i_dof,3) + 1;
      else
        i_dof = mod(i_dof-1,3) + 1;
        j_dof = 4;
      end

      dNdxii1 =  aInvJac(:,j_dof);
      dNdxij1 = -aInvJac(:,i_dof);

      dNdxii2 =  aInvJac(:,j_dof+n_vert);
      dNdxij2 = -aInvJac(:,i_dof+n_vert);

      dNdxii3 =  aInvJac(:,j_dof+2*n_vert);
      dNdxij3 = -aInvJac(:,i_dof+2*n_vert);

      vBase = cat( 3, ( aInvJac(:,i_dof+4).*dNdxii3 + aInvJac(:,j_dof+4).*dNdxij3 ) - ...
                      ( aInvJac(:,i_dof+8).*dNdxii2 + aInvJac(:,j_dof+8).*dNdxij2 ), ...
                      ( aInvJac(:,i_dof+8).*dNdxii1 + aInvJac(:,j_dof+8).*dNdxij1 ) - ...
                      ( aInvJac(:,i_dof)  .*dNdxii3 + aInvJac(:,j_dof)  .*dNdxij3 ), ...
                      ( aInvJac(:,i_dof)  .*dNdxii2 + aInvJac(:,j_dof)  .*dNdxij2 ) - ...
                      ( aInvJac(:,i_dof+4).*dNdxii1 + aInvJac(:,j_dof+4).*dNdxij1 ) );
    end

  otherwise
    vBase = 0;

end
