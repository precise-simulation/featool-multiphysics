function [ varargout ] = postplot( varargin )
%POSTPLOT Postprocessing and visualization function.
%
%   [ H ] = POSTPLOT( PROB, VARARGIN ) Postprocessing function to plot and
%   visualize FEATool problems. PROB is a valid finite element problem
%   struct. Accepts the following property/value pairs and optionally
%   returns an array H of plot object handles.
%
%       Property    Value/{Default}           1D 2D 3D Description
%       -----------------------------------------------------------------------------------
%       surfexpr    string expression         x  x  x  Surface plot expression
%       surfhexpr   string expression            x     Expression for surface height plot
%       isoexpr     string expression            x  x  Contour/iso-surface plot expression
%       isolev      scalar{10}|vector            x  x  Number/vector of contour/iso-levels
%       isores      {[20 20 20]}                    x  Isosurface interpolation resolution
%       isocolor    {'k'}                        x     Contour/iso-line color
%       arrowexpr   n_sdim string expressions    x  x  Cell with arrow plot expressions
%       arrowcolor  {r}                          x  x  Arrow color
%       arrowspacing vector {[10 10 (10)]}       x  x  Number of arrows in each direction
%       sliceexpr   string expression               x  Slice plot expression
%       slicex      [x-dir. middle of domain]       x  Slice x-position coordinates
%       slicey      [y-dir. middle of domain]       x  Slice y-position coordinates
%       slicez      [z-dir. middle of domain]       x  Slice z-position coordinates
%       deformexpr  n_sdim string expressions    x  x  Cell with deformation expressions
%       deformscale scalar {0.1}                 x  x  Deformation scaling (geom fraction)
%       solnum      scalar {n_sols}           x  x  x  Solution number/time to plot
%       evalstyle   exact|{average}           x  x     Shape function evaluation type
%       axis        off|{on}                  x  x  x  Show axes
%       grid        on |{off}                 x  x  x  Show grid
%       view        {[-35 20]}                      x  3D view setting
%       axequal     off|{on}                  x  x  x  Axis equal setting
%       edges       on |{off}                       x  Show internal edges for surface plot
%       boundary    off|{on}                     x  x  Show geometry boundary
%       bbox        {0.05}                    x  x  x  Size of bounding box (0=off)
%       alpha       {1}                             x  Transparency level
%       facecolor   {interp}                     x  x  Cell face style or fixed color
%       linestyle   {1/3D:-}|{2D:none}        x  x  x  Line style
%       linewidth   {1}                       x  x     Line width
%       edgecolor   {k}                          x     Cell edge line color
%       color       {k}                       x        Line color
%       marker      {none}                    x        Marker style
%       markersize  {5}                       x        Marker size
%       markercolor {k}                       x        Marker color
%       cellabels   on|{off}                  x  x  x  Show cell numbers
%       nodelabels  on|{off}                  x  x  x  Show vertex/node numbers
%       fontsize    {axes default}            x  x  x  Font size used in text
%       colormap    {jet}                     x  x  x  Colormap
%       colorbar    {1D:off}|{2/3D:on}        x  x  x  Colorbar toggle
%       location    {eastoutside}             x  x  x  Colorbar location
%       selexpr     string expression         x  x  x  Expression to select cells to plot
%       selcells    {all}                     x  x  x  Index vector to cells to plot
%       selsubd     {all}                     x  x  x  Index vector to subdomains to plot
%       title       string                    x  x  x  Plotly plot title
%       renderer    {empty string}|plotly     x  x  x  Render in MATLAB or Plotly
%       parent      {gca}                     x  x  x  Plot axes handle
%       setaxes     {on}|off                  x  x  x  Use axes modification settings

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help postplot, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'postplot', varargin{:} );
if( ~nargout ), clear varargout; end
