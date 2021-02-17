function [ varargout ] = intbdr( varargin )
%INTBDR Integation of expression over boundaries.
%
%   [ VAL ] = INTBDR( S_EXPR, PROB, IND_B, I_CUB, SOLNUM, IND_S )
%   Integrates the expression S_EXPR over the boundaries indicated in
%   IND_B. PROB is a valid finite element problem struct, and I_CUB
%   specifies the numerical integration rule. IND_S optionally
%   specifies which subdomain to use as reference for internal/
%   interior boundaries (normals point out from these subdomains).
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 Expression to integrate
%       prob        struct                 Finite element problem struct
%       ind_b       [1,n_bdr]              Boundary numbers (default all)
%       i_cub       scalar                 Numerical integration rule (default 2)
%       solnum      scalar {n_sols}        Solution number/time to evaluate
%       ind_s       integer array          Integration subdomains for internal boundaries
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       val         scalar                 Result of integration
%
%   See also INTSUBD, MINMAXSUBD, MINMAXBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help intbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'intbdr', varargin{:} );
if( ~nargout ), clear varargout; end
