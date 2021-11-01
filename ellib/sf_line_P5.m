function [ vBase, nLDof, xLDof, sfun ] = sf_line_P5( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_LINE_P5 1D Fifth order Lagrange shape functions for lines (P5).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LINE_P5( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming fifth order P5 Lagrange shape functions on 1D line elements
%   with values defined in the nodes and center. XI are Barycentric coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar: 1              Number of space dimensions
%       n_vert      scalar: 2              Number of vertices per cell
%       i_dof       scalar: 1-6            Local basis function to evaluate
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


nLDof = [2 0 0 4];
xLDof = [1 0 4/5 3/5 2/5 1/5;
         0 1 1/5 2/5 3/5 4/5];
sfun  = 'sf_line_P5';


switch i_eval    % Evaluation type flag.

  case 1   % Evaluation of function values.

    switch i_dof   % Basis function to evaluate.

      case 1
        vBase = (625*xi(1)^5)/24 - (625*xi(1)^4)/12 + (875*xi(1)^3)/24 - (125*xi(1)^2)/12 + xi(1);
      case 2
        vBase = - (625*xi(1)^5)/24 + (625*xi(1)^4)/8 - (2125*xi(1)^3)/24 + (375*xi(1)^2)/8 - (137*xi(1))/12 + 1;
      case 3
        vBase = - (3125*xi(1)^5)/24 + (6875*xi(1)^4)/24 - (5125*xi(1)^3)/24 + (1525*xi(1)^2)/24 - (25*xi(1))/4;
      case 4
        vBase = (3125*xi(1)^5)/12 - 625*xi(1)^4 + (6125*xi(1)^3)/12 - (325*xi(1)^2)/2 + (50*xi(1))/3;
      case 5
        vBase = - (3125*xi(1)^5)/12 + (8125*xi(1)^4)/12 - (7375*xi(1)^3)/12 + (2675*xi(1)^2)/12 - 25*xi(1);
      case 6
        vBase = (3125*xi(1)^5)/24 - (4375*xi(1)^4)/12 + (8875*xi(1)^3)/24 - (1925*xi(1)^2)/12 + 25*xi(1);
    end

  case 2   % Evaluation of first derivative.

    switch i_dof   % Basis function derivative to evaluate.

      case 1
        dNdxi1 = (3125*xi(1)^4)/24 - (625*xi(1)^3)/3 + (875*xi(1)^2)/8 - (125*xi(1))/6 + 1;
      case 2
        dNdxi1 = - (3125*xi(1)^4)/24 + (625*xi(1)^3)/2 - (2125*xi(1)^2)/8 + (375*xi(1))/4 - 137/12;
      case 3
        dNdxi1 = - (15625*xi(1)^4)/24 + (6875*xi(1)^3)/6 - (5125*xi(1)^2)/8 + (1525*xi(1))/12 - 25/4;
      case 4
        dNdxi1 = (15625*xi(1)^4)/12 - 2500*xi(1)^3 + (6125*xi(1)^2)/4 - 325*xi(1) + 50/3;
      case 5
        dNdxi1 = - (15625*xi(1)^4)/12 + (8125*xi(1)^3)/3 - (7375*xi(1)^2)/4 + (2675*xi(1))/6 - 25;
      case 6
        dNdxi1 = (15625*xi(1)^4)/24 - (4375*xi(1)^3)/3 + (8875*xi(1)^2)/8 - (1925*xi(1))/6 + 25;
    end

    vBase = aInvJac(:,1) * dNdxi1;

  case 22   % Evaluation of second derivatives.

    switch i_dof   % Basis function derivative to evaluate.

      case 1
        dNdxi1 = (125*(5*xi(1) - 2)*(10*xi(1)^2 - 8*xi(1) + 1))/12;
      case 2
        dNdxi1 = -(125*(5*xi(1) - 3)*(10*xi(1)^2 - 12*xi(1) + 3))/12;
      case 3
        dNdxi1 = - (15625*xi(1)^3)/6 + (6875*xi(1)^2)/2 - (5125*xi(1))/4 + 1525/12;
      case 4
        dNdxi1 = (15625*xi(1)^3)/3 - 7500*xi(1)^2 + (6125*xi(1))/2 - 325;
      case 5
        dNdxi1 = - (15625*xi(1)^3)/3 + 8125*xi(1)^2 - (7375*xi(1))/2 + 2675/6;
      case 6
        dNdxi1 = (15625*xi(1)^3)/6 - 4375*xi(1)^2 + (8875*xi(1))/4 - 1925/6;
    end

    vBase = -aInvJac(:,1) ./ aInvJac(:,3) * dNdxi1;

  otherwise
    vBase = 0;

end
