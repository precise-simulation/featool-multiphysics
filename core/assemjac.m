function [ varargout ] = assemjac( varargin )
%ASSEMJAC Assemble monolithic Jacobian matrix.
%
%   [ J, T_J, T_SP ] = ASSEMJAC( PROB, U, METHOD, SOLCOMP, DPERT, I_CUB, F_SPARSE )
%
%   Called from assemblenjac to assemble a monolithic Jacobian matrix J.
%
%   The METHOD flag indicates computing the Jacobian J = d(K(U)*U-f(U))/dU with respect to
%   the solution and linearization point U (METHOD=0), where each column J_i is computed as
%
%       J_i = ( (K(U_i+DPERT)*(U_i+DPERT)-f(U_i+DPERT)) - (K(U)*U-f(U)) )/DPERT
%
%   or alternatively a coordinate direction (METHOD=j), Jj = d(K(U,xj)*U-f(U,xj))/dxj where
%
%       J_i = ( (K(U,xj+DPERT)*(U,xj+DPERT)-f(U,xj+DPERT)) - (K(U,x)*U-f(U,xj)) )/DPERT
%
%   The optional input argument SOLCOMP can for (METHOD=0) be specified as a 2 x n (cell)
%   array indicating which row and column (dependent variable) combinations should be
%   computed and returned. The first column corresponds to the row and second column for
%   the blocks. SOLCOMP can also equivalently be prescribed as a cell array of pairwise
%   dependent variable name strings. For (METHOD=1/2/3) SOLCOMP can be a single vector
%   to specify the blocks/dependent variable to compute and return. If SOLCOMP is not
%   prescribed or empty all blocks and dependent variable will be computed.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Problem definition struct
%       u           vector [n_u,1]         Solution vector/linearization point, size(n_u,1)
%       method      scalar {0/i_sdim}      Jacobian computation type
%                                              0 - Jacobian with respect to solution
%                                          1/2/3 - Jacobian with respect to direction
%       solcomp     {all dvars/subd}       Dependent variables/subdomains to assemble for
%       dpert       scalar/{1e-8}          Perturbation parameter
%       icub        scalar/{2}             Numerical integration rule
%       f_sparse    scalar/{true}          Return sparse matrix format
%                                                                                         .
%       Output      Value                  Description
%       -----------------------------------------------------------------------------------
%       J           sparse/triplet         Assembled Jacobian matrix (size(n_u,n_u))
%       t_j         scalar                 Time spent assembling Jacobian matrix
%       t_sp        scalar                 Time for sparse matrix conversion
%
%   See also ASSEMBLENJAC, ASSEMBLEJ

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemjac, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemjac', varargin{:} );
if( ~nargout ), clear varargout; end
