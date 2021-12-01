function [ varargout ] = plotgeom( varargin )
%PLOTGEOM Visualization of geometry objects.
%
%   [ H ] = PLOTGEOM( SIN, VARARGIN ) Function to plot and visualize
%   geometry objects. SIN is a valid fea problem struct or cell array
%   of geometry objects. Accepts the following property/value pairs.
%
%       Property    Value/{Default}        1D 2D 3D  Description
%       -----------------------------------------------------------------------------------
%       facecolor   [.8 .8 1]              x  x  x   Face style or fixed color
%       edgecolor   {k}                    x  x  x   Edge line color
%       linestyle   {-}                    x  x  x   Edge line style
%       linewidth   {1}                    x  x  x   Edge line width
%       highlight   {[]}                   x  x  x   Highlight selected geometry objects
%       labels      off|{on}               x  x  x   Show geometry labels
%       fontsize    {2*axes default}       x  x  x   Font size used in labels
%       axequal     off|{on}               x  x  x   Axis equal setting
%       bbox        {0.05}                 x  x  x   Size of bounding box (0=off)
%       ngrid       {25}                   x  x  x   Background grid resolution
%       style       {1}                       x      Plot objects (=1) or boundaries (=2)
%       alpha       {0.5}                        x   Transparency level
%       renderer    {empty string}|plotly  x  x  x   Render in MATLAB or Plotly
%       plotlyfile  {default}              x  x  x   Filename to use for plotly export
%       parent      {gca}                  x  x  x   Plot axes handle
%
%   See also PLOTGRID, PLOTBDR, PLOTSUBD, POSTPLOT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help plotgeom, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'plotgeom', varargin{:} );
if( ~nargout ), clear varargout; end
