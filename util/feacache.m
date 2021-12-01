function [ varargout ] = feacache( varargin )
%FEACACHE Clear or examine cached expressions.
%
%   [ ST_CACHE ] = FEACACHE( CLEAR ) Clear cache (with any input
%   agument), or return a struct ST_CACHE of cached expressions used
%   by the evalexpr0, parseexpr, subcoef, and symvarbe functions.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help feacache, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'feacache', varargin{:} );
if( ~nargout ), clear varargout; end
