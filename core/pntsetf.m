function [ varargout ] = pntsetf( varargin )
%PNTSETF Assemble point source terms.
%
%   [ F, INDROW, T ] = PNTSETF( PROB, I_CUB, F )
%   Assembles point sources to the right hand side load vector F
%   with the information in the finite element problem struct PROB.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Finite element problem struct
%       i_cub       scalar                 Numerical integration rule
%       f           (neq,1)                Right hand side/load vector
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (neq,1)                Modified right hand side/load vector
%       indrow      (neq,1)                Index to rows (dofs) in rhs which were set
%       t           scalar                 Time spent in function

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help pntsetf, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'pntsetf', varargin{:} );
if( ~nargout ), clear varargout; end
