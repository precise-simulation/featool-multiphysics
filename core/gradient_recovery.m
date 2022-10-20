%GRADIENT_RECOVERY Gradient recovery with L2 projection
%
%   [ U_XJ ] = GRADIENT_RECOVERY( PROB, I_DVAR, J_DERIV, U, IND_C, I_CUB, I_MASS, XI, AINVJAC )
%
%   Computes the gradient J_DERIV (2 = x, 3 = y, 4 = z) for dependent
%   variable I_DVAR in the degrees of freedom. U is an optional
%   solution vector (default PROB.SOL.U(:,end)). IND_C specifies cells
%   to evaluate (default all), I_CUB the numerical quadrature
%   rule. I_MASS prescribes mass matrix lumping 1 = full (default for
%   higher order elements), 2 = row sum (default for linear elements),
%   3 = diagonal, 4 = HRZ lumping. When given, the optional arguments
%   XI and AINVJAC computes the derivative in local coordinates XI for
%   cells in IND_C.
%
%   See also EVALEXPR0

% Copyright 2013-2022 Precise Simulation, Ltd.
