function [ vBase, nLDof, xLDof, sfun ] = sf_lagmp( i_eval, n_sdim, n_vert, i_dof, xi, aInvJac, vBase )
%SF_LAGMP Constant lagrange multiplier.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_LAGMP( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates constant lagrange multiplier.
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

% Copyright 2013-2021 Precise Simulation, Ltd.


nLDof = [0 0 0 0];
xLDof = [];
sfun  = 'sf_lagmp';

switch i_eval   % Evaluation type flag.

  case 1        % Evaluation of function values.
    vBase = 1;

  otherwise
    vBase = 0;
end
