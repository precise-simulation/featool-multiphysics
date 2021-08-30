function [ varargout ] = postrevolve( varargin )
%POSTREVOLVE Revolve axisymmetric 2D fea data to 3D.
%
%   [ FEA ] = POSTREVOLVE( FEA, N_CR, TH ) Revolves an axisymmetric 2D
%   FEA problem struct to 3D for postprocessing. N_CR is the number of
%   cells in the direction of revolution (default 32), and TH the
%   fraction of revolution 0..1 (default 1).
%
%   Examples:
%
%      1) Revolve a axisymmetric pipe flow solution.
%
%      fea = ex_navierstokes8( 'iplot', false );
%
%      feap = postrevolve( fea, 24, 0.75 );
%
%      postplot( feap, 'surfexpr', 'sqrt(u^2+w^2)', ...
%                      'colorbar', 'off', 'axis', 'off' )
%      view(70,25)
%
%      2) Revolve a axisymmetric spherical capacitor solution.
%
%      fea = ex_electrostatics2( 'iplot', false );
%
%      feap = postrevolve( fea, 24, 0.75 );
%
%      postplot( feap, 'surfexpr', 'V', 'colorbar', 'off', 'axis', 'off' )
%      view(70,25)
%
%   See also POSTPLOT, GRIDREVOVLE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help postrevolve, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'postrevolve', varargin{:} );
if( ~nargout ), clear varargout; end
