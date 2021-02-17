function [ varargout ] = assemblea( varargin )
%ASSEMBLEA Finite element assembly of bilinear forms.
%
%   [ VROWINDS, VROWCOLS, VAVALS, N_ROWS, N_COLS ] =
%               ASSEMBLEA( FORM, C_SFUN, C_COEF, ICUB, P, C, A, SIND, N_CMAX, N_OFF, PROB, I_HRZ )
%   Assembles a matrix representation of a finite element bilinear form.
%
%   Each column in FORM specifies one bilinear additive term to compute. The
%   first row corresponds to the trial function space, and second test space.
%   The FORM entries denote either function values or derivatives to evaluate.
%
%   C_SFUN indicates which shape functions to use for the trial and test spaces.
%
%   C_COEF specifies the coefficients to evaluate for each term (column) in FORM.
%
%   The computed matrix entries are returned in VAVALS where VROWINDS and VROWCOLS
%   are corresponding row and column pointers. N_ROWS and N_COLS give the number of
%   rows in the matrix. A sparse representation of the matrix can subsequently be
%   constructed by calling SPARSE
%
%       A = sparse( vRowInds, vColInds, vAvals, n_rows, n_cols );
%
%   on the output of ASSEMLBLEA.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       form        [2,n_forms]            Forms to compute, first row indicates trial
%                                          function space (second test function space)
%                                           1     = function value
%                                           2/3/4 = x/y/z-derivatives
%       c_sfun      strings {2}            Shape functions used for the trial and
%                                          test function spaces
%       c_coef      {n_forms}              Coefficients for each bilinear form
%       i_cub       scalar                 Numerical integration rule
%       p           [n_sdim,n_p]           Array with grid point coordinates
%       c           [n_vert,n_cells]       Array with cell connectivities
%       a           [n_vert,n_cells]       Array with cell adjacency information
%       sind        [n_cells]              (optional) Subdomain indices
%       n_cmax      scalar {50000}         (optional) Max number of cells to assemble
%                                          for at once (to limit memory consumption)
%       n_off       [2]                    (optional) Row/column pointer offset
%       prob        struct                 (optional) FEATool problem struct, used
%                                          for dof mapping and evaluating coefficients
%       i_hrz       logical {false}        Apply HRZ diagonal (mass) lumping
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vRowInds    [n_A]                  Row pointers to matrix entries
%       vColInds    [n_A]                  Column pointers to matrix entries
%       vAvals      [n_A]                  Values of the assembled matrix entries
%       n_rows      scalar                 Number of rows in matrix
%       n_cols      scalar                 Number of columns in matrix
%
%   See also ASSEMBLEPROB, MAPDOFBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemblea, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemblea', varargin{:} );
if( ~nargout ), clear varargout; end
