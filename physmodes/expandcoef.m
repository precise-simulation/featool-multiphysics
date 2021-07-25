function [ varargout ] = expandcoef( varargin )
%EXPANDCOEF Expand equation/boundary coefficients.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help expandcoef, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'expandcoef', varargin{:} );
if( ~nargout ), clear varargout; end
