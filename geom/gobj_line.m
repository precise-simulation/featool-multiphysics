function [ varargout ] = gobj_line( varargin )
%GOBJ_LINE Create line geometry object.
%
%   [ GOBJ ] = GOBJ_LINE( P, TAG ) Creates a line geometry
%   object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0;1]}            Start and end line coordinates (2 x n_sdim)
%       tag         string {L1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_line, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_line', varargin{:} );
if( ~nargout ), clear varargout; end
