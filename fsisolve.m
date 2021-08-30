function [ varargout ] = fsisolve( varargin )
%FSISOLVE Solve fluid-structure interaction problem.
%
%   [ U, TLIST, P ] = FSISOLVE( PROB, VARARGIN ) Solves the fluid-
%   structure interaction problem described in the PROB finite
%   element struct. Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       tinit       scalar/{0}                   Initial time
%       tstep       scalar/{0.1}                 Time step size
%       tmax        scalar/{1}                   Maximum simulation time
%       tout        array/{tinit:tstep:tmax}     Array of output times
%       init        u0|{expr}/{0}                Initial value solution or expression
%       tnlin       string/{implicit}            Time non-linearity (semi-)implicit
%       bdf         scalar/{2}                   Time stepping scheme BDF order
%       linsolv     string/{mumps}               Sparse linear solver
%                                                    backslash, mumps, gmres, bicgstab, amg
%       maxnit      scalar/{30}                  Maximum number of non-linear iterations
%       ntol        scalar/{1e-6}                Non-linear solver tolerance
%       mmodel      string/{StVenantKirchhoff}   Material model Linear, StVenantKirchhoff,
%                                                or Neohookean (only available in 3D)
%
%   See also SOLVETIME, SOLVELIN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fsisolve, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fsisolve', varargin{:} );
if( ~nargout ), clear varargout; end
