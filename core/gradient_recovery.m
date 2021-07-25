function [ varargout ] = gradient_recovery( varargin )
%GRADIENT_RECOVERY Gradient recovery with L2 projection
%
%   [ U_XJ ] = GRADIENT_RECOVERY( PROB, I_DVAR, J_DERIV, U, IND_C, I_CUB, I_MASS, XI, AINVJAC )
%
%   Recovers the gradient J_DERIV (2-4) for dependent variable I_DVAR
%   in the degrees of freedom. U is an optional solution vector
%   (instead of taking PROB.SOL.U(:,end)). IND_C specifies cells to
%   evaluate, I_CUB the numerical quadrature rule. I_MASS prescribes
%   mass matrix lumping 1 = full (default for higher order elements),
%   2 = row sum (default for linear elements), 3 = diagonal, 4 = HRZ
%   XI and AINVJAC optionally computes the derivative in XI for cells
%   in IND_C.
%
%   See also EVALEXPR0

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gradient_recovery, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gradient_recovery', varargin{:} );
if( ~nargout ), clear varargout; end
