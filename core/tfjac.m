function [ varargout ] = tfjac( varargin )
%TFJAC Computes transformation Jacobian and determinant.
%
%   [ AJTMP, AJAC ] = TFJAC( I_FLAG, P, C, AJTMP, XI, AJAC ) Computes the inverse of
%   the transformation Jacobian and corresponding determinant in point XI
%   (in local coordinates). The input and output array AJTMP saves and reuses
%   factors that are independent of the evaluation point XI.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       i_flag      (scalar)               Evaluation flag (0=Help array, 1=All, 2=Jac.)
%       p           (n_sdim,n_p)           Array of cell vertex coordinates
%       c           (n_pcell,n_c)          Index array for cells to evaluate
%       aJTmp       (n_c,:)                Temporary help array
%       xi          (n_sdim(+1))           Local cell coordinates to evaluate
%       aJac        (n_c,:)                Array for Jacobian inverse and determinant
%       store_aJTmp (scalar)               Flag to store aJTmp as a persistent variable
%
%       Output      Value/(Size)           Description
%                                                                                         .
%       -----------------------------------------------------------------------------------
%       aJTmp       (n_c,:)                Temporary help array
%       aJac        (n_c,:)                Array for Jacobian inverse and determinant
%
%   See also TFJACLINE, TFJACTRI, TFJACQUAD, TFJACTET, TFJACHEX

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help tfjac, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'tfjac', varargin{:} );
if( ~nargout ), clear varargout; end
