function [ varargout ] = gridsmooth( varargin )
%GRIDSMOOTH Grid smoothing.
%
%   [ SOUT ] = GRIDSMOOTH( SIN, N_SM, ALPHA, I_SM ) Apply grid to
%   grid in SIN. N_SM specifies the number of smoothing steps to apply
%   (default 3) and ALPHA is the relaxation parameter (default 0.8).
%   I_SM specifies the smoothing method (default 2).
%
%   1: Laplacian smoothing operator defined as:
%
%         Pnew = (1-ALPHA)*P + ALPHA/#neighbours * SUM(Qj)
%
%   with Qj being the neighbours of the point P.
%
%   2: Density weighted umbrella operator defined as:
%
%    Pnew = (1-ALPHA) * P + ALPHA/SUM(|P-Qj|) * SUM(|P-Qj| * Qj)
%
%   with Qj being the neighbours of the point P. The constant
%   ALPHA defines the weighting of the midpoint P of the "patch".

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridsmooth, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridsmooth', varargin{:} );
if( ~nargout ), clear varargout; end
