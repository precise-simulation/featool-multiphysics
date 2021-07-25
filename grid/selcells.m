function [ varargout ] = selcells( varargin )
%SELCELLS Select cells from expression.
%
%   [ CIND ] = SELCELLS( SIN, S_EXPR ) Returns index to cells CIND
%   which fulfill the selection expression S_EXPR. If a sdim field is
%   not found in the input struct SIN the default coordinate names x,
%   y, and z will be used. MODE either indicates selection critera
%   applied to all (MODE 1, default) or any cell vertices (MODE 2).
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       sin         struct                 Grid or problem struct with
%                                          p, c (and s) fields
%       s_expr      string                 Cell selection string expression
%       mode        1/2                    Selection mode (1 = all, 2 = any vertices)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       cind        array                  Array with indicies to cells which cell
%                                          vertices satisfy the selection expression
%
%   Example:
%
%      Given a problem struct 'fea', select cells between coordinates 0.2<x<0.5:
%
%      cind = selcells( fea, '(x>0.2) & (x<0.5)' );

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help selcells, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'selcells', varargin{:} );
if( ~nargout ), clear varargout; end
