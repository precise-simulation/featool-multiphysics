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
%       dfuzz       scalar {0}                Fuzziness parameter for CSG operations
%       delint      logical {false}           Delete internal borders
%
%   The DFUZZ and DELINT parameters only apply to the default geomtool
%   engine. DFUZZ (default off) enables more tolerant (fuzzy) CSG
%   operations, collapsing close points, lines, and surfaces, so as to
%   reduce slivers and small sections. DFUZZ > 0 specifies and
%   absolute tolerance and < 0 relative to the mean geometry object
%   diameter. DELINT enables deletion of internal boundaries and
%   borders (default off).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geomcfg, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geomcfg', varargin{:} );
if( ~nargout ), clear varargout; end
