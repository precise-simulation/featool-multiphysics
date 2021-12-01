function [ varargout ] = plotgeomb( varargin )
%PLOTGEOMB Visualization of (decomposed) geometry boundaries.
%
%   [ H ] = PLOTGEOMB( SIN, VARARGIN ) Function to plot and visualize
%   (merge/decomposed) geometry object boundaries. SIN is a valid fea
%   problem struct or cell array of geometry objects. Accepts the
%   following property/value pairs.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       axis        on |{off}                 Show axes
%       grid        on |{off}                 Show grid
%       axequal     off|{on}                  Axis equal setting
%       facecolor   {[.9,.9,.9]}              2D and unselected face color
%       linestyle   {-}                       Edge line style
%       linewidth   {1}                       Edge line width
%       arrowsize   {-0.025}                  Size of arrows (negative indicates relative)
%       colors      {r,g,b,c,m,y}             Cell array with colors
%       labels      on/{off}                  Print boundary numbers
%       fontsize    {2*axes default}          Font size used in text labels
%       selbdr      {all}                     Index vector to boundaries to plot
%       plotint     {true}                    Plot internal boundaries
%       merge       {true}                    Flag to merge/decompose geometry objects
%       parent      {gca}                     Axes handle to plot in
%
%   See also PLOTGEOM, PLOTBDR, PLOTSUBD, PLOTGRID, POSTPLOT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help plotgeomb, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'plotgeomb', varargin{:} );
if( ~nargout ), clear varargout; end
