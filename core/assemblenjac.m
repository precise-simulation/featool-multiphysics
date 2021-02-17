function [ varargout ] = assemblenjac( varargin )
%ASSEMBLENJAC Assembly of defect and Newton Jacobian.
%
%   [ D, J ] = ASSEMBLENJAC( U, PROB, VARARGIN ) Evaluates defect with solution vector U
%   and optionally assembles the nonlinear Jacobian matrix.
%   Accepts optional propery value pairs in VARARGIN.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       u           vector [n_u,1]         Solution vector for linearization point
%       prob        struct                 Problem definition struct
%       jac         scalar/struct {0}      Jacobian computation type
%                                            . - Direct assembly using symbolic jac.form
%                                                and jac.coef struct contents
%                                            0 - Assembly using numerical differentiation
%                                           >0 - Directional Jacobian (jac = i_sdim dir.)
%                                           <0 - Numeric Mat-Prod Jacobian construction
%       du          scalar {1e-8}          Magnitude of numerical differentiation shift
%       icub        scalar {2}             Numerical integration rule
%       n_cmax      scalar {50000}         Max number of cells to assemble
%                                          for at once (to limit memory consumption)
%       f_sparse    scalar {1}             Return MATLAB sparse matrix format
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       d           vector [n_u,1]         Defect vector, d := f(u) - A(u)*u
%       J           sparse [n_u,n_u]       Assembled Jacobian matrix (optional output)
%
%   See also ASSEMMAT, ASSEMJAC, ASSEMBLEF

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemblenjac, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemblenjac', varargin{:} );
if( ~nargout ), clear varargout; end
