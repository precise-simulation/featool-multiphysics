function [ varargout ] = swirlflow( varargin )
%SWIRLFLOW Axisymmetric swirl flow physics mode.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help swirlflow, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'swirlflow', varargin{:} );
if( ~nargout ), clear varargout; end
