function [ varargout ] = nonnewtonian( varargin )
%NONNEWTONIAN Non-Newtonian flow physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help nonnewtonian, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'nonnewtonian', varargin{:} );
if( ~nargout ), clear varargout; end
