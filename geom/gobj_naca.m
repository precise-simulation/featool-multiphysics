function [ varargout ] = gobj_naca( varargin )
%GOBJ_NACA Create NACA 4-series wing profile geometry object.
%
%   [ GOBJ ] = GOBJ_NACA( SERIES, P0, L, TH, N, TAG ) Creates a NACA 4-series
%   wing profile geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       series      string  {0012}             NACA 4-series identifier
%       p0          array   {[0 0]}            Coordinates of leading edge/point
%       l           length  {1}                Length of wing chord
%       th          scalar  {0}                Angle of attack (degrees)
%       n           integer {20}               Number of points to define shape
%       tag         string  {N1}               Geometry object tag/name

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gobj_naca, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gobj_naca', varargin{:} );
if( ~nargout ), clear varargout; end
