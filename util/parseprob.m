function [ varargout ] = parseprob( varargin )
%PARSEPROB Check and parse a problem struct.
%
%   [ PROB ] = PARSEPROB( PROB ) Checks and parses the PROB
%   struct on the high (equation) formulation level.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Problem definition struct
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Parsed problem definition struct
%
%   See also PARSEPROB0, PARSEPHYS

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help parseprob, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'parseprob', varargin{:} );
if( ~nargout ), clear varargout; end
