function [ varargout ] = feacheckmex( varargin )
%FEACHECKMEX Check/test for valid mex file usage.
%
%   [ TF ] = FEACHECKMEX( MEX_FILE, ...ARGS ) Check/test for valid mex
%   file usage. MEX_FILE is an optional mex file to check/test against
%   (default chkmex/ftsparse), and ARGS optional input arguments in
%   test tex function call.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help feacheckmex, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'feacheckmex', varargin{:} );
if( ~nargout ), clear varargout; end
