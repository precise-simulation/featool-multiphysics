function [ vBase, nLDof, xLDof, sfun ] = sf_quad_Q3( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_QUAD_Q3 Bicubic conforming shape function for quadrilaterals (Q3).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_QUAD_Q3( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming bicubic Q3 shape functions on quadrilaterals
%   with values defined in the nodes, edges, and cell center. XI is
%   [-1..1]^2 reference coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  2             Number of space dimensions
%       n_vert      scalar:  4             Number of vertices per cell
%       i_dof       scalar: 1-16           Local basis function to evaluate
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
%   See also SFLAG3, SF_QUAD_H3

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [4 8 0 4];
xLDof = [-1  1 1 -1 -1/3  1   1/3 -1    1/3 1   -1/3 -1   -1/3  1/3 1/3 -1/3;
         -1 -1 1  1 -1   -1/3 1    1/3 -1   1/3  1   -1/3 -1/3 -1/3 1/3  1/3];
sfun  = 'sf_quad_Q3';


switch i_eval

  case 1

    switch i_dof
      case 1
        vBase = (3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) - 1);
      case 2
        vBase = -(3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) + 1)*(xi(2) - 1);
      case 3
        vBase = (3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) + 1)*(xi(2) + 1);
      case 4
        vBase = -(3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) + 1);
      case 5
        vBase = -(3*xi(1) - 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) - 1);
      case 6
        vBase = (3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((9*xi(2))/16 + 9/16)*(xi(1) + 1)*(xi(2) - 1);
      case 7
        vBase = -(3*xi(1) + 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) + 1);
      case 8
        vBase = (3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) + 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
      case 9
        vBase = (3*xi(1) + 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) - 1);
      case 10
        vBase = -(3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) + 1)*((9*xi(2))/16 + 9/16)*(xi(1) + 1)*(xi(2) - 1);
      case 11
        vBase = (3*xi(1) - 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((3*xi(2))/16 + 1/16)*(xi(1) - 1)*(xi(2) + 1);
      case 12
        vBase = -(3*xi(1) - 1)*((3*xi(1))/16 + 1/16)*(3*xi(2) - 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
      case 13
        vBase = (3*xi(1) - 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
      case 14
        vBase = -(3*xi(1) + 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) - 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
      case 15
        vBase = (3*xi(1) + 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) + 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
      case 16
        vBase = -(3*xi(1) - 1)*((9*xi(1))/16 + 9/16)*(3*xi(2) + 1)*((9*xi(2))/16 + 9/16)*(xi(1) - 1)*(xi(2) - 1);
    end

  case {2,3}

    switch i_dof
      case 1
        dNdxi1 = -((3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) - 1)*(- 27*xi(1)^2 + 18*xi(1) + 1))/256;
        dNdxi2 = -((3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) - 1)*(- 27*xi(2)^2 + 18*xi(2) + 1))/256;
      case 2
        dNdxi1 = -((3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) - 1)*(27*xi(1)^2 + 18*xi(1) - 1))/256;
        dNdxi2 = ((3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) + 1)*(- 27*xi(2)^2 + 18*xi(2) + 1))/256;
      case 3
        dNdxi1 = ((3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) + 1)*(27*xi(1)^2 + 18*xi(1) - 1))/256;
        dNdxi2 = ((3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) + 1)*(27*xi(2)^2 + 18*xi(2) - 1))/256;
      case 4
        dNdxi1 = ((3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) + 1)*(- 27*xi(1)^2 + 18*xi(1) + 1))/256;
        dNdxi2 = -((3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) - 1)*(27*xi(2)^2 + 18*xi(2) - 1))/256;
      case 5
        dNdxi1 = (9*(3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) - 1)*(- 9*xi(1)^2 + 2*xi(1) + 3))/256;
        dNdxi2 = (9*(3*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 27*xi(2)^2 + 18*xi(2) + 1))/256;
      case 6
        dNdxi1 = (9*(3*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(27*xi(1)^2 + 18*xi(1) - 1))/256;
        dNdxi2 = -(9*(3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) + 1)*(- 9*xi(2)^2 + 2*xi(2) + 3))/256;
      case 7
        dNdxi1 = -(9*(3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) + 1)*(9*xi(1)^2 + 2*xi(1) - 3))/256;
        dNdxi2 = -(9*(3*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(27*xi(2)^2 + 18*xi(2) - 1))/256;
      case 8
        dNdxi1 = -(9*(3*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 27*xi(1)^2 + 18*xi(1) + 1))/256;
        dNdxi2 = (9*(3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) - 1)*(9*xi(2)^2 + 2*xi(2) - 3))/256;
      case 9
        dNdxi1 = (9*(3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) - 1)*(9*xi(1)^2 + 2*xi(1) - 3))/256;
        dNdxi2 = -(9*(3*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(- 27*xi(2)^2 + 18*xi(2) + 1))/256;
      case 10
        dNdxi1 = -(9*(3*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(27*xi(1)^2 + 18*xi(1) - 1))/256;
        dNdxi2 = -(9*(3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) + 1)*(9*xi(2)^2 + 2*xi(2) - 3))/256;
      case 11
        dNdxi1 = -(9*(3*xi(2) - 1)*(3*xi(2) + 1)*(xi(2) + 1)*(- 9*xi(1)^2 + 2*xi(1) + 3))/256;
        dNdxi2 = (9*(3*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(27*xi(2)^2 + 18*xi(2) - 1))/256;
      case 12
        dNdxi1 = (9*(3*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(- 27*xi(1)^2 + 18*xi(1) + 1))/256;
        dNdxi2 = (9*(3*xi(1) - 1)*(3*xi(1) + 1)*(xi(1) - 1)*(- 9*xi(2)^2 + 2*xi(2) + 3))/256;
      case 13
        dNdxi1 = -(81*(3*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(- 9*xi(1)^2 + 2*xi(1) + 3))/256;
        dNdxi2 = -(81*(3*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 9*xi(2)^2 + 2*xi(2) + 3))/256;
      case 14
        dNdxi1 = -(81*(3*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(9*xi(1)^2 + 2*xi(1) - 3))/256;
        dNdxi2 = (81*(3*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(- 9*xi(2)^2 + 2*xi(2) + 3))/256;
      case 15
        dNdxi1 = (81*(3*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(9*xi(1)^2 + 2*xi(1) - 3))/256;
        dNdxi2 = (81*(3*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(9*xi(2)^2 + 2*xi(2) - 3))/256;
      case 16
        dNdxi1 = (81*(3*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 9*xi(1)^2 + 2*xi(1) + 3))/256;
        dNdxi2 = -(81*(3*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(9*xi(2)^2 + 2*xi(2) - 3))/256;
    end

    if( i_eval==2 )
      vBase = aInvJac(:,1)*dNdxi1 + aInvJac(:,2)*dNdxi2;
    else
      vBase = aInvJac(:,3)*dNdxi1 + aInvJac(:,4)*dNdxi2;
    end

  case {22,23,32,33}   % Evaluation of second order derivatives.
    error('sf_quad_Q3: second order derivative evaluation not supported.')

  otherwise
    vBase = 0;

end
