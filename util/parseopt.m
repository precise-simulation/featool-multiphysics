function [ varargout ] = parseopt( varargin )
%PARSEOPT Parse options.
%
%   [ GOT, VAL ] = PARSEOPT( COPTDEF, VARARGIN ) parses the property/value pairs
%   in VARARGIN, and adds the default options in the second column of
%   COPTDEF. The first column of COPTDEF contains strings of the field names.
%   Returns the got/val structs.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help parseopt, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'parseopt', varargin{:} );
if( ~nargout ), clear varargout; end
