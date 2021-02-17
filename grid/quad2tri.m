function [ varargout ] = quad2tri( varargin )
%QUAD2TRI Convert quadrilaterals to triangles.
%
%   [ SOUT ] = QUAD2TRI( SIN, IORDER ) Convert quadrilateral to triangular cells.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and optionally s and b) fields
%       iorder      1/2/{0}                Node ordering of triangles.
%                                            0 : Automatically select best splitting
%                                            1 : 1-2-3 / 1-3-4 splitting
%                                            2 : 2-3-4 / 1-2-4 splitting
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sout        struct                 Output grid or problem struct
%
%   See also HEX2TET, TET2HEX, TRI2QUAD

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help quad2tri, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'quad2tri', varargin{:} );
if( ~nargout ), clear varargout; end
