function [ vBase, nLDof, xLDof, sfun ] = sf_quad_Q4( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_QUAD_Q4 Biquartic conforming shape function for quadrilaterals (Q4).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_QUAD_Q4( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming biquartic Q4 shape functions on quadrilaterals
%   with values defined in the nodes, edges, and cell center. XI is
%   [-1..1]^2 reference coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  2             Number of space dimensions
%       n_vert      scalar:  4             Number of vertices per cell
%       i_dof       scalar: 1-25           Local basis function to evaluate
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
%   See also SF_QUAD_Q1

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [4 12 0 9];
xLDof = [-1  1 1 -1  -1/2  1   1/2 -1     0  1 0 -1   1/2 1   -1/2 -1    -1/2  1/2 1/2 -1/2  0   1/2 0   -1/2 0;
         -1 -1 1  1  -1   -1/2 1    1/2  -1  0 1  0  -1   1/2  1   -1/2  -1/2 -1/2 1/2  1/2 -1/2 0   1/2  0   0];
sfun  = 'sf_quad_Q4';


switch i_eval

  case 1

    switch i_dof
      case 1
        vBase = (xi(1)*xi(2)*(- 4*xi(1)^3 + 4*xi(1)^2 + xi(1) - 1)*(- 4*xi(2)^3 + 4*xi(2)^2 + xi(2) - 1))/36;
      case 2
        vBase = (xi(1)*xi(2)*(- 4*xi(1)^3 - 4*xi(1)^2 + xi(1) + 1)*(- 4*xi(2)^3 + 4*xi(2)^2 + xi(2) - 1))/36;
      case 3
        vBase = (xi(1)*xi(2)*(- 4*xi(1)^3 - 4*xi(1)^2 + xi(1) + 1)*(- 4*xi(2)^3 - 4*xi(2)^2 + xi(2) + 1))/36;
      case 4
        vBase = (xi(1)*xi(2)*(- 4*xi(1)^3 + 4*xi(1)^2 + xi(1) - 1)*(- 4*xi(2)^3 - 4*xi(2)^2 + xi(2) + 1))/36;
      case 5
        vBase = -(2*xi(1)*xi(2)*(- 2*xi(1)^3 + xi(1)^2 + 2*xi(1) - 1)*(- 4*xi(2)^3 + 4*xi(2)^2 + xi(2) - 1))/9;
      case 6
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(1)^3 - 4*xi(1)^2 + xi(1) + 1)*(- 2*xi(2)^3 + xi(2)^2 + 2*xi(2) - 1))/9;
      case 7
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(2)^3 - 4*xi(2)^2 + xi(2) + 1)*(- 2*xi(1)^3 - xi(1)^2 + 2*xi(1) + 1))/9;
      case 8
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(1)^3 + 4*xi(1)^2 + xi(1) - 1)*(- 2*xi(2)^3 - xi(2)^2 + 2*xi(2) + 1))/9;
      case 9
        vBase = -(xi(2)*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 4*xi(2)^3 + 4*xi(2)^2 + xi(2) - 1))/6;
      case 10
        vBase = -(xi(1)*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 4*xi(1)^3 - 4*xi(1)^2 + xi(1) + 1))/6;
      case 11
        vBase = -(xi(2)*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 4*xi(2)^3 - 4*xi(2)^2 + xi(2) + 1))/6;
      case 12
        vBase = -(xi(1)*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 4*xi(1)^3 + 4*xi(1)^2 + xi(1) - 1))/6;
      case 13
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(2)^3 + 4*xi(2)^2 + xi(2) - 1)*(- 2*xi(1)^3 - xi(1)^2 + 2*xi(1) + 1))/9;
      case 14
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(1)^3 - 4*xi(1)^2 + xi(1) + 1)*(- 2*xi(2)^3 - xi(2)^2 + 2*xi(2) + 1))/9;
      case 15
        vBase = -(2*xi(1)*xi(2)*(- 2*xi(1)^3 + xi(1)^2 + 2*xi(1) - 1)*(- 4*xi(2)^3 - 4*xi(2)^2 + xi(2) + 1))/9;
      case 16
        vBase = -(2*xi(1)*xi(2)*(- 4*xi(1)^3 + 4*xi(1)^2 + xi(1) - 1)*(- 2*xi(2)^3 + xi(2)^2 + 2*xi(2) - 1))/9;
      case 17
        vBase = (16*xi(1)*xi(2)*(- 2*xi(1)^3 + xi(1)^2 + 2*xi(1) - 1)*(- 2*xi(2)^3 + xi(2)^2 + 2*xi(2) - 1))/9;
      case 18
        vBase = (16*xi(1)*xi(2)*(- 2*xi(2)^3 + xi(2)^2 + 2*xi(2) - 1)*(- 2*xi(1)^3 - xi(1)^2 + 2*xi(1) + 1))/9;
      case 19
        vBase = (16*xi(1)*xi(2)*(- 2*xi(1)^3 - xi(1)^2 + 2*xi(1) + 1)*(- 2*xi(2)^3 - xi(2)^2 + 2*xi(2) + 1))/9;
      case 20
        vBase = (16*xi(1)*xi(2)*(- 2*xi(1)^3 + xi(1)^2 + 2*xi(1) - 1)*(- 2*xi(2)^3 - xi(2)^2 + 2*xi(2) + 1))/9;
      case 21
        vBase = (4*xi(2)*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 2*xi(2)^3 + xi(2)^2 + 2*xi(2) - 1))/3;
      case 22
        vBase = (4*xi(1)*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 2*xi(1)^3 - xi(1)^2 + 2*xi(1) + 1))/3;
      case 23
        vBase = (4*xi(2)*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 2*xi(2)^3 - xi(2)^2 + 2*xi(2) + 1))/3;
      case 24
        vBase = (4*xi(1)*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 2*xi(1)^3 + xi(1)^2 + 2*xi(1) - 1))/3;
      case 25
        vBase = (4*xi(1)^4 - 5*xi(1)^2 + 1)*(4*xi(2)^4 - 5*xi(2)^2 + 1);
    end

  case {2,3}

    switch i_dof
      case 1
        dNdxi1 = -(xi(2)*(4*xi(1) - 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(- 4*xi(1)^2 + 2*xi(1) + 1))/36;
        dNdxi2 = -(xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) - 1)*(xi(1) - 1)*(- 4*xi(2)^2 + 2*xi(2) + 1))/36;
      case 2
        dNdxi1 = (xi(2)*(4*xi(1) + 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(4*xi(1)^2 + 2*xi(1) - 1))/36;
        dNdxi2 = -(xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) - 1)*(xi(1) + 1)*(- 4*xi(2)^2 + 2*xi(2) + 1))/36;
      case 3
        dNdxi1 = (xi(2)*(4*xi(1) + 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) + 1)*(4*xi(1)^2 + 2*xi(1) - 1))/36;
        dNdxi2 = (xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) + 1)*(xi(1) + 1)*(4*xi(2)^2 + 2*xi(2) - 1))/36;
      case 4
        dNdxi1 = -(xi(2)*(4*xi(1) - 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) + 1)*(- 4*xi(1)^2 + 2*xi(1) + 1))/36;
        dNdxi2 = (xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) + 1)*(xi(1) - 1)*(4*xi(2)^2 + 2*xi(2) - 1))/36;
      case 5
        dNdxi1 = (2*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(- 8*xi(1)^3 + 3*xi(1)^2 + 4*xi(1) - 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) - 1)*(4*xi(2) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 4*xi(2)^2 + 2*xi(2) + 1))/9;
      case 6
        dNdxi1 = -(2*xi(2)*(4*xi(1) + 1)*(2*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(4*xi(1)^2 + 2*xi(1) - 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(xi(1) + 1)*(- 8*xi(2)^3 + 3*xi(2)^2 + 4*xi(2) - 1))/9;
      case 7
        dNdxi1 = (2*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) + 1)*(- 8*xi(1)^3 - 3*xi(1)^2 + 4*xi(1) + 1))/9;
        dNdxi2 = -(2*xi(1)*(2*xi(1) + 1)*(4*xi(2) + 1)*(xi(1) - 1)*(xi(1) + 1)*(4*xi(2)^2 + 2*xi(2) - 1))/9;
      case 8
        dNdxi1 = (2*xi(2)*(4*xi(1) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 4*xi(1)^2 + 2*xi(1) + 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(xi(1) - 1)*(- 8*xi(2)^3 - 3*xi(2)^2 + 4*xi(2) + 1))/9;
      case 9
        dNdxi1 = (xi(1)*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(8*xi(1)^2 - 5)*(xi(2) - 1))/3;
        dNdxi2 = -((2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 4*xi(2)^2 + 2*xi(2) + 1))/6;
      case 10
        dNdxi1 = ((4*xi(1) + 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(4*xi(1)^2 + 2*xi(1) - 1))/6;
        dNdxi2 = (xi(1)*xi(2)*(2*xi(1) - 1)*(2*xi(1) + 1)*(8*xi(2)^2 - 5)*(xi(1) + 1))/3;
      case 11
        dNdxi1 = (xi(1)*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(8*xi(1)^2 - 5)*(xi(2) + 1))/3;
        dNdxi2 = ((2*xi(1) - 1)*(2*xi(1) + 1)*(4*xi(2) + 1)*(xi(1) - 1)*(xi(1) + 1)*(4*xi(2)^2 + 2*xi(2) - 1))/6;
      case 12
        dNdxi1 = -((4*xi(1) - 1)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 4*xi(1)^2 + 2*xi(1) + 1))/6;
        dNdxi2 = (xi(1)*xi(2)*(2*xi(1) - 1)*(2*xi(1) + 1)*(8*xi(2)^2 - 5)*(xi(1) - 1))/3;
      case 13
        dNdxi1 = (2*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) - 1)*(- 8*xi(1)^3 - 3*xi(1)^2 + 4*xi(1) + 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) + 1)*(4*xi(2) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 4*xi(2)^2 + 2*xi(2) + 1))/9;
      case 14
        dNdxi1 = -(2*xi(2)*(4*xi(1) + 1)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(4*xi(1)^2 + 2*xi(1) - 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(xi(1) + 1)*(- 8*xi(2)^3 - 3*xi(2)^2 + 4*xi(2) + 1))/9;
      case 15
        dNdxi1 = (2*xi(2)*(2*xi(2) - 1)*(2*xi(2) + 1)*(xi(2) + 1)*(- 8*xi(1)^3 + 3*xi(1)^2 + 4*xi(1) - 1))/9;
        dNdxi2 = -(2*xi(1)*(2*xi(1) - 1)*(4*xi(2) + 1)*(xi(1) - 1)*(xi(1) + 1)*(4*xi(2)^2 + 2*xi(2) - 1))/9;
      case 16
        dNdxi1 = (2*xi(2)*(4*xi(1) - 1)*(2*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(- 4*xi(1)^2 + 2*xi(1) + 1))/9;
        dNdxi2 = (2*xi(1)*(2*xi(1) - 1)*(2*xi(1) + 1)*(xi(1) - 1)*(- 8*xi(2)^3 + 3*xi(2)^2 + 4*xi(2) - 1))/9;
      case 17
        dNdxi1 = -(16*xi(2)*(2*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(- 8*xi(1)^3 + 3*xi(1)^2 + 4*xi(1) - 1))/9;
        dNdxi2 = -(16*xi(1)*(2*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 8*xi(2)^3 + 3*xi(2)^2 + 4*xi(2) - 1))/9;
      case 18
        dNdxi1 = -(16*xi(2)*(2*xi(2) - 1)*(xi(2) - 1)*(xi(2) + 1)*(- 8*xi(1)^3 - 3*xi(1)^2 + 4*xi(1) + 1))/9;
        dNdxi2 = -(16*xi(1)*(2*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(- 8*xi(2)^3 + 3*xi(2)^2 + 4*xi(2) - 1))/9;
      case 19
        dNdxi1 = -(16*xi(2)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 8*xi(1)^3 - 3*xi(1)^2 + 4*xi(1) + 1))/9;
        dNdxi2 = -(16*xi(1)*(2*xi(1) + 1)*(xi(1) - 1)*(xi(1) + 1)*(- 8*xi(2)^3 - 3*xi(2)^2 + 4*xi(2) + 1))/9;
      case 20
        dNdxi1 = -(16*xi(2)*(2*xi(2) + 1)*(xi(2) - 1)*(xi(2) + 1)*(- 8*xi(1)^3 + 3*xi(1)^2 + 4*xi(1) - 1))/9;
        dNdxi2 = -(16*xi(1)*(2*xi(1) - 1)*(xi(1) - 1)*(xi(1) + 1)*(- 8*xi(2)^3 - 3*xi(2)^2 + 4*xi(2) + 1))/9;
      case 21
        dNdxi1 = -(8*xi(1)*xi(2)*(2*xi(2) - 1)*(8*xi(1)^2 - 5)*(xi(2) - 1)*(xi(2) + 1))/3;
        dNdxi2 = (4*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 8*xi(2)^3 + 3*xi(2)^2 + 4*xi(2) - 1))/3;
      case 22
        dNdxi1 = (4*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 8*xi(1)^3 - 3*xi(1)^2 + 4*xi(1) + 1))/3;
        dNdxi2 = -(8*xi(1)*xi(2)*(2*xi(1) + 1)*(8*xi(2)^2 - 5)*(xi(1) - 1)*(xi(1) + 1))/3;
      case 23
        dNdxi1 = -(8*xi(1)*xi(2)*(2*xi(2) + 1)*(8*xi(1)^2 - 5)*(xi(2) - 1)*(xi(2) + 1))/3;
        dNdxi2 = (4*(4*xi(1)^4 - 5*xi(1)^2 + 1)*(- 8*xi(2)^3 - 3*xi(2)^2 + 4*xi(2) + 1))/3;
      case 24
        dNdxi1 = (4*(4*xi(2)^4 - 5*xi(2)^2 + 1)*(- 8*xi(1)^3 + 3*xi(1)^2 + 4*xi(1) - 1))/3;
        dNdxi2 = -(8*xi(1)*xi(2)*(2*xi(1) - 1)*(8*xi(2)^2 - 5)*(xi(1) - 1)*(xi(1) + 1))/3;
      case 25
        dNdxi1 = 2*xi(1)*(8*xi(1)^2 - 5)*(4*xi(2)^4 - 5*xi(2)^2 + 1);
        dNdxi2 = 2*xi(2)*(8*xi(2)^2 - 5)*(4*xi(1)^4 - 5*xi(1)^2 + 1);
    end

    if( i_eval==2 )
      vBase = aInvJac(:,1).*dNdxi1 + aInvJac(:,2).*dNdxi2;
    else
      vBase = aInvJac(:,3).*dNdxi1 + aInvJac(:,4).*dNdxi2;
    end

  case {22,23,32,33}   % Evaluation of second order derivatives.
    error('sf_quad_Q4: second order derivative evaluation not supported.')

  otherwise
    vBase = 0;

end
