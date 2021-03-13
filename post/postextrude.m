function [ varargout ] = postextrude( varargin )
%POSTEXTRUDE Extrude 2D fea data to 3D.
%
%   [ FEA ] = POSTEXTRUDE( FEA, NZ, LZ, AX ) Extrudes a 2D FEA problem
%   struct to 3D for postprocessing. NZ is the number of cells in the
%   direction of extrusion (default 5), and LZ the corresponding
%   extrusion length (default 1). AX [-3..3] optionally rotates and
%   aligns the extruded grid with the selected axis (default -2 which
%   is the negative y-axis). Non-vertex ordered dependent variable
%   solutions will be interpolated to the grid point vertices.
%
%   Examples:
%
%      1) Extrude and plot a backwards facing step flow field.
%
%      fea = ex_navierstokes4( 'iplot', false );
%
%      feap = postextrude( fea );
%
%      postplot( feap, 'isoexpr','sqrt(u^2+v^2)', 'sliceexpr', 'p', ...
%                      'slicex', [], 'slicey', -sqrt(eps), 'slicez', [], ...
%                      'colorbar', 'off', 'axis', 'off' )
%
%      2) Extrude and plot a deflected beam.
%
%      fea = ex_piezoelectric1( 'iplot', false );
%
%      feap = postextrude( fea, 5, 6e-3 );
%
%      postplot( feap, 'surfexpr', 'u', 'deformexpr', {'u',0,'v'}, ...
%                      'deformscale', 2, 'colorbar', 'off', 'axis', 'off' )
%      plotgeom( feap, 'labels', 'off' )
%
%      axis normal
%      axis( [ 0, 0.121, -20e-3, 20e-3, -20e-3, 20e-3 ] )
%      view( 33, 24 )
%
%   See also POSTPLOT, GRIDEXTRUDE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help postextrude, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'postextrude', varargin{:} );
if( ~nargout ), clear varargout; end
