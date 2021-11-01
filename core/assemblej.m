function [ varargout ] = assemblej( varargin )
%ASSEMBLEJ Finite element assembly of Jacobian forms.
%
%   [ VROWINDS, VROWCOLS, VJVALS, N_ROWS, N_COLS ] = ASSEMBLEJ( AFORM, FFORM, C_SFUN,
%        C_ACOEF, C_FCOEF, ICUB, P, C, A, IND_S, PROB, OFFSET, IPERT, DPERT, N_CMAX )
%
%   Assembles a numerical Jacobian matrix representation of a finite element
%   bilinear and linear forms.
%
%   IPERT indicates to perturb the solution vector (=0) with DPERT (default 1e-8)
%   calculating the Jacobian J = dK/du*U - f, where each column J_i is given as
%
%       J_i = ( (K(U_i+DPERT)*(U_i+DPERT)-f(U_i+DPERT)) - (K(U)*U-f(U)) )/DPERT
%
%   Alternatively, if IPERT is equal to an integer (1, 2 or 3) the Jacobian with
%   respect to space dimension/coordinate IPERT will be computed as
%
%       J_i = ( (K(U,x+DPERT)*(U,x+DPERT)-f(U,x+DPERT)) - (K(U,x)*U-f(U,x)) )/DPERT
%
%   The linearization point/solution vector U needs to be given and stored in
%   the prob.sol.u field. If it is not present a zero solution will be used.
%
%   Each column in AFORM and FFORM specify bilinear and linear additive terms to
%   compute. In AFORM the first row corresponds to the trial function space, and
%   second test space. The form entries denote either function values or
%   derivatives to evaluate.
%
%   C_SFUN indicates which shape functions to use for the trial and test function
%   spaces. Optionally, a third entry is required for the trial function space when
%   assembling off-diagonal (target) blocks for coupled problems. If only one or two
%   C_SFUN entries is provided, the first entry will be used for target trial
%   function space.
%
%   For off-diagonal blocks the dof OFFSET ( dof_offset_row, dof_offset_column )
%   is used to find and perturb the target block solution vector.
%
%   C_ACOEF and C_FCOEF specifies the coefficients to evaluate for the forms.
%
%   The computed matrix entries are returned in VJVALS where VROWINDS and VROWCOLS
%   are corresponding row and column pointers. N_ROWS and N_COLS give the number of
%   rows in the matrix. J sparse representation of the matrix can subsequently be
%   constructed by calling SPARSE
%
%       J = sparse( vRowInds, vColInds, vJvals, n_rows, n_cols );
%
%   on the output of ASSEMLBLEJ.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       aform       [2,n_aforms]           Bilinear forms to compute, first row indicates
%                                          trial function space (second test function space)
%                                           1     = function value
%                                           2/3/4 = x/y/z -derivatives
%       fform       [n_fforms]             Linear forms to compute
%       c_sfun      strings {2/3}          Shape functions used for the assmebly trial,
%                                          assembly test, and target trial function spaces
%       c_acoef     {n_aforms}             Coefficients for each bilinear form
%       c_fcoef     {n_fforms}             Coefficients for each linear form
%       i_cub       scalar                 Numerical integration rule
%       p           [n_sdim,n_p]           Array with grid point coordinates
%       c           [n_vert,n_cells]       Array with cell connectivities
%       a           [n_vert,n_cells]       Array with cell adjacency information
%       ind_s       [n_cells]              (optional) Subdomain indices
%       prob        struct                 (optional) FEATool problem struct, used
%       offset      [2]                    (optional) Offset for finding start in solution
%                                          for matrix coef/vector product dof perturbation
%       ipert       scalar/{0/i_sdim}      0-Solution perturbation, >0 Coordinate dir. pert
%       dpert       scalar                 (optional) Perturbation parameter (default 1e-8)
%       n_cmax      scalar {50000}         Max number of cells per group to assemble
%                                          for at once (to limit memory consumption)
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       vRowInds    [n_J]                  Row pointers to matrix entries
%       vColInds    [n_J]                  Column pointers to matrix entries
%       vJvals      [n_J]                  Values of the assembled matrix entries
%       n_rows      scalar                 Number of rows in matrix
%       n_cols      scalar                 Number of columns in matrix
%
%   See also ASSEMJAC, MAPDOFBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemblej, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemblej', varargin{:} );
if( ~nargout ), clear varargout; end
