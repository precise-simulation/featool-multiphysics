function [ varargout ] = intsubd( varargin )
%INTSUBD Integate expression over subdomains.
%
%   [ VAL ] = INTSUBD( S_EXPR, PROB, IND_S, IND_C, I_CUB, SOLNUM ) Integrates
%   the expression S_EXPR over the subdomains indicated in IND_S, or
%   alternatively the cells in IND_C. PROB is a valid finite element
%   problem struct, and I_CUB specifies the numerical integration rule.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 Expression to integrate
%       prob        struct                 Finite element problem struct
%       ind_s       [1,n_subd]             Subdomain numbers (default all)
%       ind_c       [1,n_cells]            Cell indices (default all)
%       i_cub       scalar                 Numerical integration rule (default 2)
%       solnum      scalar {n_sols}        Solution number/time to evaluate
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       val         scalar                 Result of integration
%
%   See also INTBDR, MINMAXSUBD, MINMAXBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help intsubd, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'intsubd', varargin{:} );
if( ~nargout ), clear varargout; end
