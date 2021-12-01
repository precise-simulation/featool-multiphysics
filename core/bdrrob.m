function [ varargout ] = bdrrob( varargin )
%BDRROB Compute Robin boundary condition contributions.
%
%   [ VROWINDS, VROWCOLS, VAVALS, N_ROWS, N_COLS ] =
%               BDRROB( PROB, I_CUB ) Computes matrix contributinos
%   for Robin boundary conditions to the problem specified in the
%   BDR.N field of the problem struct PROB. It is assumbed the problem
%   has been linearized by usin BDR_LIN.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Finite element problem struct
%       i_cub       scalar                 Numerical integration rule
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       vRowInds    [n_A]                  Row pointers to matrix entries
%       vColInds    [n_A]                  Column pointers to matrix entries
%       vAvals      [n_A]                  Values of the assembled matrix entries
%       n_rows      scalar                 Number of rows in matrix
%       n_cols      scalar                 Number of columns in matrix
%
%   See also BDRNEU_LIN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help bdrrob, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'bdrrob', varargin{:} );
if( ~nargout ), clear varargout; end
