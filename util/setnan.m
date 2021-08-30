function [ varargout ] = setnan( varargin )
%SETNAN Set nan values.
%
%   [ X ] = SETNAN( X, VAL ) Sets all nan values in the numeric array
%   X to VAL (default 0).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help setnan, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'setnan', varargin{:} );
if( ~nargout ), clear varargout; end
