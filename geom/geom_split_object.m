function [ varargout ] = geom_split_object( varargin )
%GEOM_SPLIT_OBJECT Split composite geometry object.

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom_split_object, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom_split_object', varargin{:} );
if( ~nargout ), clear varargout; end
