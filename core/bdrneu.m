%BDRNEU Compute Neumann (flux) boundary condition vector.
%
%   [ F, T ] = BDRNEU( PROB, I_CUB ) Computes a load vector F with
%   Neumann/flux boundary condition contributions for the boundary
%   expressions defined in PROB.BDR.N.
%
%   Boundary coefficients are specified as (1 x n_dvar) cell arrays in
%   the BDR.N field where the entry for each boundary is a (n_bc_groups
%   x n_bdr) nested cell array containing the coefficients
%   (n_bc_groups is 1, except for special element types such as
%   Hermite basis functions). The coefficient entries can be specified
%   either as constant numeric values, or string expressions which
%   will be evaluated during the simulation.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Problem data struct
%       i_cub       scalar {2}             Numerical integration rule
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (neq,1)                Neumann boundary condition load vector
%       t           scalar                 Time spent in function

% Copyright 2013-2026 Precise Simulation, Ltd.
