function [ varargout ] = fenics( varargin )
%FENICS Exports, solves, and imports a FEniCS project.
%
%   [ PROB ] = FENICS( PROB, VARARGIN ) Export, solves, or imports the solved
%   problem described in the PROB finite element struct using the FEniCS project
%   solver. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       modes       data, export, solve, import  Command mode(s) to call (default all)
%       data        default                      FEniCS input data file to use
%       fname       featool-fenics               FEniCS base filename root
%       fdir        string  {pwd}                Directory to write mesh and data files
%       order       scalar  {2}                  Integration order used in c-expressions
%       scmd        system default               Custom system solve command string
%       ccmd        string                       FEniCS installation bash check command
%                   python -c "import dolfin;print(dolfin.dolfin_version())"
%       wbash       C:\Windows\System32\bash.exe Windows bash executable path
%       clear       boolean {true}               Clear output and log files
%       pause       boolean {false}              Pause external log terminal
%       fid         scalar  {1}                  File identifier for output ([]=no output)
%
%   Examples:
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

% Copyright 2013-2019 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fenics, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fenics', varargin{:} );
if( ~nargout ), clear varargout; end
