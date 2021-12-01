function [ vBase, nLDof, xLDof, sfun ] = sf_line_H3( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_LINE_H3 Third order 1D C1 Hermite shape functions for lines.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_H3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates C1 Hermite shape functions on 1D line elements with value and
%   first derivatives defined in the nodes. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-4            Local basis function to evaluate
%       xi          array  [2,1]           Local coordinates of evaluation point
%       aInvJac     [n,3]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [2,4]                  Number of local degrees of freedom on
%                                          vertices, edges, faces, cell interiors,
%                                          and vertices without boundary conditions
%       xLDof       [2,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_LINE_P3

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [2 0 0 0;
         2 0 0 0];
xLDof = [1 0 1 0;
         0 1 0 1];
sfun  = 'sf_line_H3';


switch i_eval    % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = 3*xi(1)^2 - 2*xi(1)^3;
      case 2
        vBase = 3*xi(2)^2 - 2*xi(2)^3;
      case 3
        vBase = ( xi(1)^2 - xi(1)^3 ) .* aInvJac(:,3);
      case 4
        vBase = ( xi(2)^3 - xi(2)^2 ) .* aInvJac(:,3);
    end


  case 2   % Evaluation of first derivatives.

    switch i_dof   % Basis function derivative to evaluate.
      case 1
        dNdxi1 = -6*xi(1)*(xi(1) - 1);
        dNdxi2 =  0;
      case 2
        dNdxi1 =  0;
        dNdxi2 = -6*xi(2)*(xi(2) - 1);
      case 3
        dNdxi1 = -xi(1)*(3*xi(1) - 2) .* aInvJac(:,3);
        dNdxi2 =  0;
      case 4
        dNdxi1 =  0;
        dNdxi2 =  xi(2)*(3*xi(2) - 2) .* aInvJac(:,3);
    end

    vBase = aInvJac(:,1) .* dNdxi1 + aInvJac(:,2) .* dNdxi2;


  case 22   % Evaluation of second derivatives.

    switch i_dof   % Basis function derivative to evaluate.
      case 1
        dNdxi1 =  6*(2*xi(1)-1) ./ aInvJac(:,3);
        dNdxi2 =  0;
      case 2
        dNdxi1 =  0;
        dNdxi2 = -6*(2*xi(2)-1) ./ aInvJac(:,3);
      case 3
        dNdxi1 =  6*xi(1) - 2;
        dNdxi2 =  0;
      case 4
        dNdxi1 =  0;
        dNdxi2 =  6*xi(2) - 2;
    end

    vBase = aInvJac(:,1) .* dNdxi1 + aInvJac(:,2) .* dNdxi2;


  otherwise
    vBase = 0;

end
