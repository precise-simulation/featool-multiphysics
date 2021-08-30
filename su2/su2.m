function [ varargout ] = su2( varargin )
%SU2 MATLAB SU2 CFD solver CLI interface.
%
%   [ U, TLIST, VARS ] = SU2( PROB, VARARGIN ) Export, solves, or
%   imports the solved problem described in the PROB finite element
%   struct using the SU2 CFD solver. Accepts the following
%   property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       mode        check, export, solve, import Command mode(s) to call (default all)
%       cfg         default                      Default SU2 config
%       turb        scalar {0}                   Turbulence model:      0: none/laminar
%                                                  1: Spalart-Allmaras, 2:k-Omega (SST)
%       init        scalar {[]}                  Initial values   []: init expressions
%                                                  i: use solution
%       nproc       scalar  {numcores/2}         Number of processors to use
%       workdir     default                      SU2 work directory
%       fname       default                      SU2 work filename
%       logfname    default                      SU2 log/output filename
%       fid/logfid  scalar {1}                   Log file/message output file handle
%       hax         handle                       Axis handle to plot convergence
%       naxts       scalar {250}                 Maximum number of time steps to plot
%
%   Also accepts the following SU2 cfg property/value pairs to set the
%   cfg file during export
%
%       Property      Value/{Default}         Description
%       -----------------------------------------------------------------------------------
%       fname         default                 SU2 work filename
%       solver        string {INC_NAVIER_STOKES}   SU2 governing equations/solver
%       ischeme       scalar/{0}              Time stepping scheme
%                                                  0 - Stationary
%                                                  1 - Time stepping
%                                                  2 - Dual time stepping (1st order)
%                                                  3 - Dual time stepping (2nd order)
%       tstep         scalar/{0.1}            Time step size
%       tmax          scalar/{1}              Maximum simulation time
%       wrtfreq       scalar/{1}              Frequency to output solution files/to screen
%       tol           scalar {1e-8}           Stopping criteria for (P) residuals
%       maxit         scalar/{9999/20}        Maximum number of iterations (stat/timedep)
%       rho           scalar {1.0}            Density (constant)
%       miu           scalar {1.0}            Viscosity (constant)
%       upwind        string {venk_wang/muscl} Discretization scheme, valid options are
%                                             the central schemes JST and LAX-FRIEDRICH,
%                                             or UPWIND (equivalent to FDS for incompressible,
%                                             ROE for compressible) with MUSCL and slope limiting
%                                             NONE, VENKATAKRISHNAN, VENKATAKRISHNAN_WANG,
%                                             BARTH_JESPERSEN, or VAN_ALBADA_EDGE
%       init          vector {0,0,0}          Initial values
%       restart       string {}               Restart (CSV) file
%       mesh          string {mesh.su2}       SU2 mesh filename
%       nsdim         int {2}                 Number of space dimensions
%       isaxi         boolean {false}         Axisymmetric mode (only compressible)
%
%   Examples:
%
%      1) Laminar Hagen-Poiseuille flow in a channel with convergence plot.
%
%      n = 20; rho = 1; miu = 1; uin = 1;
%
%      fea.sdim = {'x' 'y'};
%      fea.geom.objects = { gobj_rectangle(0,3,0,1) };
%      fea.grid = rectgrid( 3*n, 1*n, [0 3;0 1] );
%
%      fea = addphys(fea,@navierstokes);
%      fea.phys.ns.eqn.coef{1,end} = { rho };
%      fea.phys.ns.eqn.coef{2,end} = { miu };
%      fea.phys.ns.eqn.coef{5,end} = { uin };
%      fea.phys.ns.bdr.sel(2) = 4;
%      fea.phys.ns.bdr.sel(4) = 2;
%      fea.phys.ns.bdr.coef{2,end}{1,4} = uin;
%
%      fea = parsephys( fea );
%      fea = parseprob( fea );
%
%      fea.sol.u = su2( fea, 'tol', 1e-6, 'hax', axes() );
%
%      figure
%      subplot(2,1,1)
%      postplot( fea, 'surfexpr', 'p', 'isoexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'} )
%
%      subplot(2,1,2), hold on, grid on
%      xlabel('Velocity profile at outlet'), ylabel('y')
%      x = 3*ones(1,100);
%      y = linspace(0,1,100);
%      U_ref = 6*uin*(y.*(1-y))./1^2;
%      U = evalexpr( 'sqrt(u^2+v^2)', [x;y], fea );
%      plot( U_ref, y, 'r--', 'linewidth', 3 )
%      plot( U, y, 'b-', 'linewidth', 2.5 )
%      legend( 'Analytic solution', 'Computed solution' )
%
%   See also EX_NAVIERSTOKES1 -7,10,12, EX_COMPRESSIBLEEULER2 -4

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help su2, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'su2', varargin{:} );
if( ~nargout ), clear varargout; end
