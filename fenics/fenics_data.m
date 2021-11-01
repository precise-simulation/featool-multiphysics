function [ varargout ] = fenics_data( varargin )
%FENICS_DATA FEATool to FEniCS data conversion.
%
%   [ S ] = FENICS_DATA( PROB, VARARGIN ) Converts the FEATool FEA
%   problem struct PROB to a FEniCS Python simulation script. If the
%   output S is requested the simulation script is returned as a
%   string, otherwise it is written to the file designated by the FDIR
%   and FNAME arguments. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       fname       featool-fenics               FEniCS base filename root
%       fdir        string  {}                   Directory to write mesh and data files
%       order       scalar  {2}                  Integration order used in c-expressions
%       linsolv     string  {default}            Linear solver
%       precond     string  {default}            Preconditioner for iterative solver
%       ischeme     scalar  {0}                  Solver type/time stepping scheme
%                                                  0 - Steady state/stationary
%                                                  1 - 1st order Backward Euler
%                                                  2 - 2nd order Crank-Nicolson
%       tstep       scalar  {0.1}                Time step size (for ischeme >= 1)
%       tmax        scalar  {1.0}                Maximum simulation time (for ischeme >= 1)
%       islin       scalar  {auto}               Specify linear or non-linear solver
%       maxnit      scalar  {20}                 Maximum number of non-linear iterations
%       nlrlx       scalar  {1.0}                Relaxation for non-linear iterations
%       toldef      scalar  {1e-6}               Relative stopping criteria
%       tolchg      scalar  {1e-6}               Absolute stopping criteria
%       saveall     logical {true}               Save all time steps/solutions,
%                                                otherwise just saves last one
%   Example:
%
%      1) Convert the ex_navierstokes1 example to a FEniCS script.
%
%      fea = ex_navierstokes1();   % Set up, solve and return fea data struct.
%      s = fenics_data( fea )      % Return equivalent FEniCS script as string.
%
%   See also FENICS

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fenics_data, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fenics_data', varargin{:} );
if( ~nargout ), clear varargout; end
