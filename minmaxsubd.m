function [ varargout ] = minmaxsubd( varargin )
%MINMAXSUBD Calculate subdomain minima and maxima.
%
%   [ MIN_VAL, MAX_VAL, MIN_COORD, MAX_COORD ] = MINMAXSUBD( S_EXPR, PROB, IND_S, IND_C, I_CUB, SOLNUM )
%   Evaluates the mimimum and maximum value of expression S_EXPR over the
%   subdomains indicated in IND_S or alternatively the cells in IND_C.
%   PROB is a valid finite element problem struct. Returns the minima and
%   maxima in MIN_VAL and MAX_VAL, and the corresponding coordinates in
%   MIN_COORD and MAX_COORD.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 Expression to evaluate
%       prob        struct                 Finite element problem struct
%       ind_s       [1,n_subd]             Subdomain numbers (default all)
%       ind_c       [1,n_cells]            Cell indices (default all)
%       i_cub       scalar                 Evaluation point rule (default 2)
%       solnum      scalar {n_sols}        Solution number/time to evaluate
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       min_val     scalar                 Minimum value of expression
%       max_val     scalar                 Maximum value of expression
%
%   See also MINMAXBDR, INTSUBD, INTBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help minmaxsubd, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'minmaxsubd', varargin{:} );
if( ~nargout ), clear varargout; end
