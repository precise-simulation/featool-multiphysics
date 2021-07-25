function [ varargout ] = minmaxbdr( varargin )
%MINMAXBDR Calculate boundary minima and maxima.
%
%   [ MIN_VAL, MAX_VAL, MIN_COORD, MAX_COORD ] = MINMAXBDR( S_EXPR, PROB, IND_B, I_CUB, SOLNUM )
%   Evaluates the mimimum and maximum value of expression S_EXPR over
%   the boundaries indicated in IND_B. PROB is a valid finite element
%   problem struct. Returns the minima and maxima in MIN_VAL and MAX_VAL,
%   and the corresponding coordinates in MIN_COORD and MAX_COORD.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       s_expr      string                 Expression to evaluate
%       prob        struct                 Finite element problem struct
%       ind_b       [1,n_bdr]              Boundary numbers (default all)
%       i_cub       scalar                 Evaluation point rule (default 2)
%       solnum      scalar {n_sols}        Solution number/time to evaluate
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       min_val     scalar                 Minimum value of expression
%       max_val     scalar                 Maximum value of expression
%
%   See also MINMAXSUBD, INTSUBD, INTBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help minmaxbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'minmaxbdr', varargin{:} );
if( ~nargout ), clear varargout; end
