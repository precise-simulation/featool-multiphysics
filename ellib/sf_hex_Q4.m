%SF_HEX_Q4 Triquartic conforming shape function for hexahedrons (Q4).
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = SF_HEX_Q4( I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates conforming triquartic Q4 shape functions on hexahedrons with values defined
%   in the nodes, edges, faces, and cell center. XI is [-1..1]^3 reference coordinates.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       i_eval      scalar:  1             Evaluate function values
%                           >1             Evaluate values of derivatives
%       n_sdim      scalar:  3             Number of space dimensions
%       n_vert      scalar:  8             Number of vertices per cell
%       i_dof       scalar: 1-n_ldof       Local basis function to evaluate
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
%   See also SFLAG4, SF_HEX_Q1

% Copyright 2013-2023 Precise Simulation, Ltd.

















