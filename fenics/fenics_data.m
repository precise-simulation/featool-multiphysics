function [ varargout ] = fenics_data( varargin )
%FENICS_DATA FEATool to FEniCS data conversion.
%
%   [ S ] = FENICS_DATA( PROB, VARARGIN ) Converts the FEATool
%   FEA problem struct PROB to a FEniCS Python simulation script.
%   If output S is requested the simulation script is returned as a
%   string, otherwise it is written to the file designated by
%   the FDIR and FNAME property/value pairs. Accepts the following
%   property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       fname       featool-fenics               FEniCS base filename root
%       fdir        string  {}                   Directory to write mesh and data files
%       order       scalar  {2}                  Integration order used in c-expressions
%       ischeme     scalar  {0}                  Solver type/time stepping scheme
%                                                  0 - Steady state/stationary
%                                                  1 - 1st order Backward Euler
%                                                  2 - 2nd order Crank-Nicolson
%       tstep       scalar  {0.1}                Time step size
%       tmax        scalar  {1.0}                Maximum simulation time
%       maxnit      scalar  {20}                 Maximum number of non-linear iterations
%       nlrlx       scalar  {1.0}                Relaxation for non-linear iterations
%       toldef      scalar  {1e-6}               Relative stopping criteria
%       tolchg      scalar  {1e-6}               Absolute stopping criteria
%       saveall     logical {true}               Save all time steps/solutions,
%                                                otherwise just saves last one
%   Example:
%
%      1) Convert ex_navierstokes1 example to a FEniCS script.
%
%      fea = ex_navierstokes1();
%      s = fenics_data( fea )
%
%   See also FENICS, IMPEXP_DOLFIN

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fenics_data, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fenics_data', varargin{:} );
if( ~nargout ), clear varargout; end
