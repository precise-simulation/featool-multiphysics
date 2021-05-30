function [ varargout ] = tri2quad( varargin )
%TRI2QUAD Convert triangles to quadrilaterals.
%
%   [ SOUT ] = TRI2QUAD( SIN ) Convert triangular to quadrilateral cells.
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
%   See also HEX2TET, QUAD2TRI, TET2HEX

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help tri2quad, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'tri2quad', varargin{:} );
if( ~nargout ), clear varargout; end
