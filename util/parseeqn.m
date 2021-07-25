function [ varargout ] = parseeqn( varargin )
%PARSEEQN Parse string equations and return eqn struct.
%
%   [ EQN ] = PARSEEQN( CEQN, CDVAR, CSDIM ) Parse string equation(s) and return
%   weak finite element equation struct form and coefficients.
%
%       Input        Value/Type             Description
%       -----------------------------------------------------------------------------------
%       ceqn         cell array/string      Equation string or cell array of strings
%       cdvar        cell array/string      Dependent variable name(s)
%       csdim        cell array/string      Space dimension coordinate name(s)
%
%   String quations supports the following syntax
%
%       u`     - time derivative for u (` = single quote)
%       u      - dependent variable u (explicit)
%       u_t    - dependent variable u (implicit)
%       ux     - space derivative for u in x-direction (explicit, rhs)
%       u_x    - space derivative for u in x-direction (implicit)
%       ux_t   - space derivative for u in x-direction (implicit)
%       ux_x   - 2nd space derivative for u in x-direction
%       (.)_t  - Multiplication of expression with test function
%       (.)_x  - Multiplication of expr. with test function in x-direction
%       (.)_tx - Multiplication of expr. with test function in x-direction
%
%
%   an underscore indicates multiplication with a test function, t, operating
%   on the preceding expression. For expressions in parentheses with multiple
%   dependent variable factors the operation is applied to the first from the
%   reversed direction, for example (2*uy+u*ux)_x will be expanded to
%   2*uy_x+u*ux_x. Application of derivatives _x to constants, variables, or
%   first of dependent variables will also change the equation sign accounting
%   partial integration in the weak fem formulations. In addition expressions
%   can be built up with standard MATLAB commands, functions, and operators
%   such + - * / sqrt() log() etc. For elements supporting the curl operator a
%   corresponding _c postfix is also supported.
%
%   Example:
%
%      1) Poisson's equation in 1D:
%
%      seqn = 'da*u'' - (D*ux)_x = f';
%      eqn  = parseeqn( seqn, 'u', 'x' );
%
%      Gives, eqn.m: the mass matrix, in this case scaled by da
%
%      eqn.m.form{:} >> [1;1]
%      eqn.m.coef{:} >> da
%
%      eqn.a: implicit bilinear form, here diffusion matrix scaled by D
%
%      eqn.a.form{:} >> [2;2]
%      eqn.a.coef{:} >> D
%
%      eqn.f: linear right hand side source term
%
%      eqn.f.form{:} >> 1
%      eqn.f.coef{:} >> f
%
%   See also EXPANDEQN, ASSEMBLEA, ASSEMBLEF

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help parseeqn, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'parseeqn', varargin{:} );
if( ~nargout ), clear varargout; end
