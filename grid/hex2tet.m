function [ varargout ] = hex2tet( varargin )
%HEX2TET Convert hexahedra to tetrahedra.
%
%   [ SOUT ] = HEX2TET( SIN, MODE ) Convert hexahedral to tetrahedral cells.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and optionally s and b) fields
%       mode        scalar                  1 = split each hex cell into 12 tets (default)
%                                           2 = split each hex cell into 6 tets
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Output grid or problem struct
%
%   See also QUAD2TRI, TET2HEX, TRI2QUAD

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help hex2tet, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'hex2tet', varargin{:} );
if( ~nargout ), clear varargout; end
