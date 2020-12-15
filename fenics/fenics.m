function [ varargout ] = fenics( varargin )
%FENICS MATLAB FEniCS project FEA solver CLI interface.
%
%   [ PROB ] = FENICS( PROB, VARARGIN ) Export, solves, or imports the solved
%   problem described in the PROB finite element struct using the FEniCS project
%   solver. The fea problem struct is returned with PTOB.SOL.U and PROB.SOL.T
%   fields if applicable. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       modes       {c, e, s, i, l}              Command mode(s) to call (default all)
%                                             (c)heck, (e)xport, (s)olve, (i)mport, c(l)ear
%       data        default                      FEniCS input data file to use
%       fname       featool-fenics               FEniCS base filename root
%       fdir        string  {pwd}                Directory to write mesh and data files
%       ask         boolean {false}              Ask for file name for export
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
%       saveall     logical {true}               Save all time steps/solutions
%                                                otherwise just saves last one
%       scmd        system default               Custom system solve command string
%       ccmd        string  {dolfin-version}     FEniCS installation bash check command
%       wbash       C:\Windows\System32\bash.exe Windows bash executable path
%       python      string  {auto}               Python interpreter call command
%       pause       boolean {false}              Pause external log terminal
%       fid         scalar  {1}                  File identifier for output ([]=no output)
%
%   Example:
%
%      1) Uniform heating of a unit circle.
%
%      fea.sdim = {'x', 'y'};
%      fea.grid = quad2tri(circgrid());
%
%      fea = addphys(fea, @heattransfer);
%      fea.phys.ht.eqn.coef{6,end} = {1};   % Heat source term.
%      fea.phys.ht.bdr.sel(:) = 1;          % Zero temperature boundary conditions.
%
%      fea = parsephys(fea);                % Parse physics mode and fea problem struct.
%      fea = parseprob(fea);
%
%      fea = fenics(fea);
%
%      postplot(fea, 'surfexpr', 'T')
%
%   See also FENICS_DATA, IMPEXP_DOLFIN, EX_HEATTRANSFER1-9, EX_NAVIERSTOKES1-4, EX_POISSON8

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fenics, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fenics', varargin{:} );
if( ~nargout ), clear varargout; end
