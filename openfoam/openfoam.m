function [ varargout ] = openfoam( varargin )
%OPENFOAM MATLAB OpenFOAM CFD solver CLI interface.
%
%   [ U, TLIST, VARS ] = OPENFOAM( PROB, VARARGIN ) Export, solves,
%   and/or imports the solved problem described in the PROB finite
%   element struct using the OpenFOAM CFD solver. Accepts the
%   following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       mode        check, export, solve, import Command mode(s) to call (default all)
%       data        default                      Default OpenFOAM data and parameter dict
%       init        scalar {[]}                  Initial values   []: init expressions
%                                                  i: use solution -1: potential flow
%       turb        struct                       Turbulence data struct fields
%                                                  model, inlet [k,e/o], wallfcn [1/0]
%       interp      scalar {2}                   Interpolate solution to grid points
%                                                  1: no weighting 2: cell volume weighting
%       control     logical {false}              Show solver control panel.
%       casedir     default                      OpenFOAM case directory
%       foamdir     default                      OpenFOAM installation directory
%       logfname    default                      OpenFOAM log/output filename
%       fid/logfid  scalar {1}                   Log file/message output file handle
%       hax         handle                       Axis handle to plot convergence
%       pmaxts      integer {500}                Maximum number of time steps to print/plot
%
%   MODE is a string or cell array of strings selecting action(s) to
%   perform. By default check, export, solve, and import are performed
%   in sequence.
%
%   INIT by default takes the initial value expressions from a
%   Navier-Stokes or Compressible Euler physics mode. It can also be
%   overridden to select a previous solution (integer i specifies the
%   solution number), or use potentialFoam to compute the initial
%   values.
%
%   TURB is a struct with fields, turb.model indicating turbulence
%   model (kEpsilon, realizableKE RNGkEpsilon, kOmega, kOmegaSST, or
%   SpalartAllmaras), turb.inlet a vector with two components
%   specifying the inlet values k/epsilon or k/omega when using the
%   corresponding models (can also be computed using the
%   turbulence_inletbccalc function), and turb.wallfcn is a logical
%   flag designating if turbulent wall functions should be used.
%
%   Also accepts the following OpenFOAM data property/value pairs to
%   set the control dicts during export
%
%       Property       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       application    string {simpleFoam}    OpenFOAM application binary to run
%       ddtScheme      string {steadyState}   Time stepping scheme
%       tolres         scalar/vector {1e-4}   Stopping criteria for residuals (Simple)
%       startTime      scalar {0.0}           Simulation start time
%       endTime        scalar {1000}          Simulation end time
%       deltaT         scalar {1.0}           Time step size
%       maxDeltaT      scalar {0.1}           Maximum time step size
%       maxCo          scalar {0.5}           Maximum Courant number
%       upwind         string {linearUpwind}  Discretization scheme (div), valid options are
%                                             linear, LUST, linearUpwind, limitedLinear, and upwind
%       bound          scalar {2}             Bounding/limiting, >0 bound div, >1 bound grad
%       ortho          scalar {auto}          Correction for grid non-orthogonality
%       nproc          scalar {numcores/2}    Number of processors to use
%       writeInterval  scalar {endTime}       Solution output write interval
%       writePrecision scalar {6}             Output file write precision
%       timePrecision  scalar {6}             Time format precision
%       purgeWrite     scalar {0}             Specified number of output solutions
%       transportModel string {Newtonian}     Transport model
%       simulationType string {laminar}       Simulation type laminar/RAS/LES
%       RASModel       string {kEpsilon}      RAS turbulence model: kEpsilon, realizableKE
%                                          RNGkEpsilon, kOmega, kOmegaSST, SpalartAllmaras
%       nu             scalar {1.0}           Kinematic viscosity (constant)
%
%   Examples:
%
%      1) Laminar Hagen-Poiseuille flow in a channel.
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
%      fea.sol.u = openfoam( fea );
%
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
%      2) Axisymmetric turbulent flow in a pipe, showing solution convergence curves.
%
%      Re = 1e5; rho = 1; miu = 1/Re; win = 1;
%
%      fea.sdim = {'r' 'z'};
%      fea.geom.objects = { gobj_rectangle(0,0.5,0,15) };
%      fea.grid = rectgrid( 0.5-[0 0.01 0.03 0.06 0.1 0.3 0.5], 50, [0 0.5;0 15] );
%      fea.grid = gridrefine( fea.grid );
%
%      fea = addphys(fea,{@navierstokes,true});
%      fea.phys.ns.eqn.coef{1,end} = { rho };
%      fea.phys.ns.eqn.coef{2,end} = { miu };
%      fea.phys.ns.eqn.coef{6,end} = { win };
%      fea.phys.ns.bdr.sel(1) = 2;
%      fea.phys.ns.bdr.sel(2) = 1;
%      fea.phys.ns.bdr.sel(3) = 4;
%      fea.phys.ns.bdr.sel(4) = 5;
%      fea.phys.ns.bdr.coef{2,end}{2,1} = win;
%
%      fea = parsephys( fea );
%      fea = parseprob( fea );
%
%      turb.model = 'kEpsilon';
%      turb.inlet = [0.001,0.00045];
%      turb.wallfcn = 1;
%      fea.sol.u = openfoam( fea, 'hax', axes(), 'control', true, 'turb', turb );
%
%      figure,subplot(2,1,1)
%      postplot( fea, 'surfexpr', 'p', 'isoexpr', 'sqrt(u^2+w^2)', 'arrowexpr', {'u' 'w'} )
%      axis([0,0.5,14,15])
%
%      subplot(2,1,2), hold on, grid on
%      xlabel('Velocity profile at outlet'), ylabel('r')
%      r = linspace(0,0.5,100);
%      z = 15*ones(1,100);
%      U = evalexpr( 'sqrt(u^2+w^2)', [r;z], fea );
%      plot( U, r, 'b-', 'linewidth', 2.5 )
%
%   See also TURBULENCE_INLETBCCALC, EX_NAVIERSTOKES1 -8,10-13, EX_COMPRESSIBLEEULER2 -4

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help openfoam, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'openfoam', varargin{:} );
if( ~nargout ), clear varargout; end
