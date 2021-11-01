function [ varargout ] = bdrsetd( varargin )
%BDRSETD Set Dirichlet boundary conditions.
%
%   [ AMAT, F, T, INDDOF ] = BDRSETD( AMAT, F, PROB, ISYMM, SET_NULL, SOLCOMP )
%   Sets Dirichlet boundary conditions in the global matrix AMAT and
%   right hand side/load vector F with the information in the finite
%   element problem struct PROB. The optional flag ISYMM symmetrizes
%   the boundary conditions. SET_NULL prescribes zeros instead of the
%   expression specified in the prob.bdr.d field. INDDOF is an index
%   array to the degrees of freedom corresponding to Dirichlet constraints.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       amat        (n_a,n_a)              System matrix (sparse or triplet format)
%       f           (neq,1)                Right hand side/load vector
%       prob        struct                 Finite element problem struct
%       isymm       scalar/{0}             Symmetrize BCs if applicable
%       set_null    scalar/{0}             Set zeros in f vector
%       solcomp     {all dvars/subd}       Dependent variables/subdomains to set BCs for
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       amat        (n_a,n_a)              Modified system matrix
%       f           (neq,1)                Modified right hand side/load vector
%       t           scalar                 Time spent in function
%       inddof      array                  Dof index array indices to Dirichlet constraints

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help bdrsetd, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'bdrsetd', varargin{:} );
if( ~nargout ), clear varargout; end
