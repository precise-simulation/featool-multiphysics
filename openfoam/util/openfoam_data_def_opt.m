%OPENFOAM_DATA_DEF_OPT Default OpenFOAM data options.
%
%   [ COPTDEF ] = OPENFOAM_DATA_DEF_OPT() Default OpenFOAM data options.
%
%       Property       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       application    string {simpleFoam}    OpenFOAM solver application to run
%                                               simpleFoam, pimpleFoam, rhoCentralFoam,
%                                               chtMultiRegion(Simple)Foam
%       units          string {SI}            Unit system (SI or USCS)
%       ddtSchemes     string {steadyState}   Time stepping scheme
%                                               steady state: steadyState, euler, localEuler
%                                               time dependent: backward, CrankNicolson (0.9)
%       tolres         scalar {1e-4/U, 1e-2/p}/vector   Stopping criteria for residuals
%       startTime      scalar {0.0}           Simulation start time
%       endTime        scalar {1000.0}        Simulation end time
%       deltaT         scalar {1.0}           Time step size
%       adjustTimeStep string no/{yes}        Adjust time step according to max Courant number in transient simulation
%       maxCo          scalar {0.2}           Maximum Courant number
%       maxDi          scalar {10.0}          Maximum Diffusion number (only for cht)
%       maxDeltaT      scalar {1.0}           Maximum time step size
%       upwind         string {linearUpwind}  Discretization scheme (div), valid options are
%                                               linear, LUST, linearUpwind, limitedLinear, upwind
%       bound          scalar {2}             Bounding/limiting, >0 bound div, >1 bound grad
%       ortho          scalar {0.33}          Correction for grid non-orthogonality
%       writeControl   string {adjustableRunTime}  Controls the timing of write output to file
%       writeInterval  scalar {endTime}       Solution output write interval
%       purgeWrite     scalar {0/1 steady}    Specified number of output solutions
%       writePrecision scalar {8}             Output file write precision
%       timePrecision  scalar {6}             Time format precision
%       transportModel string {Newtonian}     Transport model
%       simulationType string {laminar}       Simulation type laminar, RAS, or LES
%       RASModel       string {kEpsilon}      RAS turbulence model kEpsilon, realizableKE
%                                          RNGkEpsilon, kOmega, kOmegaSST, SpalartAllmaras
%       LESModel       string {Smagorinsky}   LES turbulence model
%       turb           struct                 Turbulence data struct fields
%                                                  model, inlet [k,e/o], wallfcn [1/0]
%       init           scalar {[]}            Initial values   []~ use init expressions
%                                                  i~ use solution i -1~ potential flow
%
%   TURB is a struct with fields, turb.model indicating turbulence
%   model (kEpsilon, realizableKE RNGkEpsilon, kOmega, kOmegaSST, or
%   SpalartAllmaras), turb.inlet a vector with two components
%   specifying the inlet values k/epsilon or k/omega when using the
%   corresponding models (can also be computed using the
%   turbulence_inletbccalc function), and turb.wallfcn is a logical
%   flag designating if turbulent wall functions should be used.
%
%   INIT by default takes the initial value expressions from the
%   physics modes. It can also be overridden to select a previous
%   solution (integer init = i specifies the solution number), or use
%   potentialFoam to compute the initial values (init = -1).
%
%   See also OPENFOAM_DATA

% Copyright 2013-2026 Precise Simulation, Ltd.
