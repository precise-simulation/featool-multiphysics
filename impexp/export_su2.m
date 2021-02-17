function [ varargout ] = export_su2( varargin )
%EXPORT_SU2 Export grid in SU2 format.
%
%   [ GRID ] = EXPORT_SU2( FILE_NAME, DATA ) Export of SU2 grid
%   data. FILE_NAME is a string specifying the (root) file name to
%   process. DATA can be either a full fea struct or just the grid
%   data struct.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help export_su2, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'export_su2', varargin{:} );
if( ~nargout ), clear varargout; end
