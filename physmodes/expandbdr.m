function [ varargout ] = expandbdr( varargin )
%EXPANDBDR Expand physics mode boundary coefficients.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help expandbdr, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'expandbdr', varargin{:} );
if( ~nargout ), clear varargout; end
