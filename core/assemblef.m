function [ varargout ] = assemblef( varargin )
%ASSEMBLEF Finite element assembly of linear form (rhs/load vector).
%
%   [ F ] = ASSEMBLEF( FORM, C_SFUN, C_COEF, ICUB, P, C, A, SIND, N_CMAX, PROB )
%   Assembles a finite element linear form (right hand side/load vector).
%   Each column in FORM specify one linear additive term to compute. The FORM
%   entries denote either function values or derivatives to evaluate. C_COEF
%   specifies the coefficients to evaluate for each term in FORM. C_SFUN
%   indicates which shape functions to use for the evaluation.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       form        [n_forms]              Linear forms to compute and add
%                    1                      - function value
%                    2,3,4                  - x,y,z -derivatives
%       c_sfun      string                 Shape functions used for the trial space
%       c_coef      [n_forms]              Coefficients to multiply/scale each form
%       i_cub       scalar                 Numerical integration rule
%       p           [n_sdim,n_p]           Array with grid point coordinates
%       c           [n_vert,n_cells]       Array with cell connectivities, for each
%                                          cell (column) points to vertices in p
%       a           [n_vert,n_cells]       Array with cell adjacency info, for each
%                                          cell (column) points to cell neighbours
%       sind        [n_cells]              Subdomain indices
%       n_cmax      scalar {50000}         (optional) Max number of cells to assemble
%                                          for at once (to limit memory consumption)
%       prob        struct                 Finite element problem struct, (optionally)
%                                          used for evaluating the coefficients
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       f    [n_f]                         Right hand side/load vector
%
%   See also ASSEMBLEPROB, MAPDOFBDR

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemblef, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemblef', varargin{:} );
if( ~nargout ), clear varargout; end
