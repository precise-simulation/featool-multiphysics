function [ varargout ] = tet2hex( varargin )
%TET2HEX Convert tetrahedra to hexahedra.
%
%   [ SOUT ] = TET2HEX( SIN ) Convert tetrahedral to hexahedral cells.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and optionally s and b) fields
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Output grid or problem struct
%
%   See also HEX2TET, QUAD2TRI, TRI2QUAD

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help tet2hex, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'tet2hex', varargin{:} );
if( ~nargout ), clear varargout; end
