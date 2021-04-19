function [ varargout ] = assembleprob( varargin )
%ASSEMBLEPROB Assembly of system matrix and right hand side/load vector.
%
%   [ M, A, F, T_M, T_A, T_F, T_SP ] = ASSEMBLEPROB( PROB, VARARGIN ) Calls the
%   assembly routines to compute and assemble a monolithic sparse system
%   matrix and a right hand side/load vector as specified in the PROB struct.
%   Accepts optional propery value pairs in VARARGIN.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Problem definition struct
%       icub        scalar/{auto}          Numnerical integration rule used in assembly
%                                                Default 1+max(shape function order)
%       imass       scalar {1}             Mass matrix lumping:  1 = Full mass matrix
%                                          2 = row sum lumping,  3 = diagonal lumping
%                                          4 = HRZ diagonal lumping
%       n_cmax      scalar {50000}         Max number of cells to assemble
%                                          for at once (to limit memory consumption)
%       f_m         scalar {0}             Assembly flag for mass matrix
%       f_a         scalar {0}             Assembly flag for system matrix
%       f_f         scalar {0}             Assembly flag for load vector
%       f_c         scalar {1}             Assembly flag for integral constraints
%       f_sparse    scalar {0}             Return MATLAB sparse matrix format
%       solcomp     {all dvars/subd}       Dependent variables/subdomains to assemble for
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       M           sparse [n_M]           Assembled mass matrix
%       A           sparse [n_A]           Assembled system matrix
%       f           [n_A,1]                Assembled rhs/load vector
%       t_m         scalar                 Time spent assembling mass matrix
%       t_a         scalar                 Time spent assembling system matrix
%       t_f         scalar                 Time spent assembling rhs/load vector
%       t_sp        scalar                 Time for sparse matrix conversion
%
%   See also ASSEMBLEA, ASSEMBLEF, ASSEMMAT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assembleprob, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assembleprob', varargin{:} );
if( ~nargout ), clear varargout; end
