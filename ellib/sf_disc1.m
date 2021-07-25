function [ vBase, nLDof, xLDof, sfun ] = sf_disc1( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_DISC1 Linear discontinuous shape function (P1).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_DISC1( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates discontinuous linear shape functions with with values
%   defined by a constant and derivatives in the cell interior.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1-3            Number of space dimensions
%       n_vert      scalar: 2-8            Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
%       xi          [n_sdim(+1)]           Local coordinates of evaluation point
%       aInvJac     [n,n_sdim(+1)*n_sdim]  Inverse of transformation Jacobian
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
%   See also SF_DISC0

% Copyright 2013-2021 Precise Simulation, Ltd.


% Check for evaluation on simplicies.
if( n_vert==n_sdim+1 )
  f_simplex = 1;
else
  f_simplex = 0;
end


nLDof = [0 0 0 n_sdim+1];
if( f_simplex )
  xLDof = ones(n_sdim+1,1)/(n_sdim+1);
else
  xLDof = zeros(n_sdim,1);
end
xLDof = repmat(xLDof,1,n_sdim+1);
sfun  = 'sf_disc1';


switch i_eval   % Evaluation type flag.

  case {1}   % Evaluation of function values.

    if ( f_simplex || i_dof<=n_sdim )

      vBase = xi(i_dof);

    else

      vBase = 1;
    end

  case {2,3,4}   % Evaluation of first order derivatives.

    if ( f_simplex || i_dof<=n_sdim )

      i_col = i_dof+(i_eval-2)*(n_sdim+f_simplex);

      vBase = aInvJac(:,i_col);

    else

      vBase = 0;
    end

  otherwise
    vBase = 0;

end
