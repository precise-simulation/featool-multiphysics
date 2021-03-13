function [ varargout ] = turbulent_mixing_length( varargin )
%TURBULENT_MIXING_LENGTH Computes turbulent mixing length.
%
%   [ L_MIX ] = TURBULENT_MIXING_LENGTH( FEA, KAPPA, C, WALL_BDR, VARARGIN )
%
%   Computes the turbulent mixing length for algebraic turbulence models
%
%       L_MIX = min( KAPPA*L_WALL, C*L_MAX )
%
%   where KAPPA (default 0.41) and C (default 0.09) are constants,
%   L_WALL the shortest distance from an evaluation point to the walls
%   indicated in WALL_BDR, and L_MAX the maximum distance of L_WALL.
%   The evaluation coordinates are given as a cell array of 2-3 column
%   vectors in VARARGIN. Computed results will be cached and reused.

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help turbulent_mixing_length, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'turbulent_mixing_length', varargin{:} );
if( ~nargout ), clear varargout; end
