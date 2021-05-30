function [ vBase, nLDof, xLDof, sfun ] = sf_line_P3( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_LINE_P3 1D Third order Lagrange shape functions for lines (P3).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_P3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming third order P3 Lagrange shape functions on 1D line elements
%   with values defined in the nodes and center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-4            Local basis function to evaluate
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
%   See also SFLAG3, SF_LINE_H3

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [2 0 0 2];
xLDof = [1 0 2/3 1/3;
         0 1 1/3 2/3];
sfun  = 'sf_line_P3';


switch i_eval    % Evaluation type flag.

  case 1   % Evaluation of function values.
    xi = xi(1);

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = xi*(2 - 3*xi)*(1 - 3*xi)/2;
      case 2
        vBase = (1 - xi)*(2 - 3*xi)*(1 - 3*xi)/2;
      case 3
        vBase = 9*xi*(1-xi)*(3*xi - 1)/2;
      case 4
        vBase = 9*xi*(1-xi)*(2 - 3*xi)/2;
    end


  case 2   % Evaluation of first derivative.
    xi = xi(1);

    switch i_dof   % Basis function derivative to evaluate.
      case 1
        dNdxi = (27*xi^2)/2 - 9*xi + 1;
      case 2
        dNdxi = 18*xi - (27*xi^2)/2 - 11/2;
      case 3
        dNdxi = 36*xi - (81*xi^2)/2 - 9/2;
      case 4
        dNdxi = (81*xi^2)/2 - 45*xi + 9;
    end

    vBase = aInvJac(:,1)*dNdxi;

  case 22   % Evaluation of second derivatives.
    xi = xi(1);

    switch i_dof   % Basis function derivative to evaluate.
      case 1
        dNdxi = 27*xi - 9;
      case 2
        dNdxi = 18 - 27*xi;
      case 3
        dNdxi = 36 - 81*xi;
      case 4
        dNdxi = 81*xi - 45;
    end

    vBase = -aInvJac(:,1)./aInvJac(:,3)*dNdxi;

  otherwise
    vBase = 0;

end
