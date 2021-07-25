function [ varargout ] = princse( varargin )
%PRINCSE Compute principal stresses and strains.
%
%   [ SE_P ] = PRINCSE( SE_X, SE_Y, SE_Z, SE_XY, SE_YZ, SE_XZ, ICOMP ) Computes
%   principal stresses and strains (for the ICOMP component, default 1). In two
%   dimensions only inputs SE_X, SE_Y, SE_Z, SE_XY, and ICOMP are required.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help princse, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'princse', varargin{:} );
if( ~nargout ), clear varargout; end
