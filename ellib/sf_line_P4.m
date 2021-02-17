function [ vBase, nLDof, xLDof, sfun ] = sf_line_P4( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_LINE_P4 1D Fourth order Lagrange shape functions for lines (P4).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_P4( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fourth order P4 Lagrange shape functions on 1D line elements
%   with values defined in the nodes and center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-5            Local basis function to evaluate
%       xi          array [2,1]            Local coordinates of evaluation point
%       aInvJac     [n,3]                  Inverse of transformation Jacobian
%       vBase       [n]                    Preallocated output vector
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vBase       [n]                    Evaluated function values
%       nLDof       [4]                    Number of local degrees of freedom on
%                                          vertices, edges, faces, and cell interiors
%       xLDof       [2,n_ldof]             Local coordinates of local dofs
%       sfun        string                 Function name of called shape function
%
%   See also SF_LINE_P1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [2 0 0 3];
xLDof = [1 0 3/4 1/2 1/4;
         0 1 1/4 1/2 3/4];
sfun  = 'sf_line_P4';


switch i_eval    % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = (32*xi(1)^4)/3 - 16*xi(1)^3 + (22*xi(1)^2)/3 - xi(1);
      case 2
        vBase = (32*xi(1)^4)/3 - (80*xi(1)^3)/3 + (70*xi(1)^2)/3 - (25*xi(1))/3 + 1;
      case 3
        vBase = - (128*xi(1)^4)/3 + (224*xi(1)^3)/3 - (112*xi(1)^2)/3 + (16*xi(1))/3;
      case 4
        vBase = 64*xi(1)^4 - 128*xi(1)^3 + 76*xi(1)^2 - 12*xi(1);
      case 5
        vBase = - (128*xi(1)^4)/3 + 96*xi(1)^3 - (208*xi(1)^2)/3 + 16*xi(1);
    end

  case 2   % Evaluation of first derivative.

    switch i_dof   % Basis function derivative to evaluate.

      case 1
        dNdxi1 = ((8*xi(1) - 3)*(16*xi(1)^2 - 12*xi(1) + 1))/3;
      case 2
        dNdxi1 = ((8*xi(1) - 5)*(16*xi(1)^2 - 20*xi(1) + 5))/3;
      case 3
        dNdxi1 = - (512*xi(1)^3)/3 + 224*xi(1)^2 - (224*xi(1))/3 + 16/3;
      case 4
        dNdxi1 = 4*(2*xi(1) - 1)*(32*xi(1)^2 - 32*xi(1) + 3);
      case 5
        dNdxi1 = - (512*xi(1)^3)/3 + 288*xi(1)^2 - (416*xi(1))/3 + 16;
    end

    vBase = aInvJac(:,1) * dNdxi1;

  case 22   % Evaluation of second derivatives.

    switch i_dof   % Basis function derivative to evaluate.

      case 1
        dNdxi1 = 128*xi(1)^2 - 96*xi(1) + 44/3;
      case 2
        dNdxi1 = 128*xi(1)^2 - 160*xi(1) + 140/3;
      case 3
        dNdxi1 = - 512*xi(1)^2 + 448*xi(1) - 224/3;
      case 4
        dNdxi1 = 768*xi(1)^2 - 768*xi(1) + 152;
      case 5
        dNdxi1 = - 512*xi(1)^2 + 576*xi(1) - 416/3;
    end

    vBase = -aInvJac(:,1) ./ aInvJac(:,3) * dNdxi1;

  otherwise
    vBase = 0;

end
