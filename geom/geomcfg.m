function [ varargout ] = geomcfg( varargin )
%GEOMCFG Sets and retrieves geometry configuration parameters.
%
%   [ VAL ] = GEOMCFG( S_NAME, VAL ) Loads the geom parameters and
%   sets up the CFG struct. S_NAME is the config target field name and
%   the optional argument VAL is passed to set the value. Call without
%   input arguments returns the whole config struct. Accepts the
%   following parameters.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       engine2     string {built-in}         2D engine (built-in, geomtool)
%       binary2     string {}                 Location of 2D geometry engine binary
%       engine3     string {built-in}         3D engine (built-in, geomtool, gmsh, brl)
%       binary3     string {}                 Location of 3D geometry engine binary

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geomcfg, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geomcfg', varargin{:} );
if( ~nargout ), clear varargout; end
