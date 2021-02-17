function [ varargout ] = bdrneu( varargin )
%BDRNEU Compute Neumann (flux) boundary condition vector.
%
%   [ F, T ] = BDRNEU( PROB, I_CUB ) Computes a (flux) load vector with Neumann boundary
%   condition contributions to the problem specified in the problem struct PROB.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Finite element problem struct
%       i_cub       scalar                 Numerical integration rule
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (neq,1)                Neumann boundary condition load vector
%       t           scalar                 Time spent in function

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help bdrneu, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'bdrneu', varargin{:} );
if( ~nargout ), clear varargout; end
