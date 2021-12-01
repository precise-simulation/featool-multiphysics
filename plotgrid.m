function [ varargout ] = plotgrid( varargin )
%PLOTGRID Grid plot function.
%
%   [ H ] = PLOTGRID( PROB, VARARGIN ) Function to plot and visualize
%   a grid. PROB is a valid grid or problem struct.
%   Accepts the following property/value pairs.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       axis        off|{on}                  Show axes
%       grid        on |{off}                 Show grid
%       axequal     off|{on}                  Axis equal setting
%       edges       on |{off}                 Show internal edges (3D only)
%       view        {[45 45]}                 3D view setting
%       bbox        {0.05}                    Size of bounding box (0=off)
%       alpha       {1}                       Transparency level
%       facecolor   {[.9 .9 .9]}              Cell face style or fixed color
%       linestyle   {-}                       Line style
%       linewidth   {1}                       Line width
%       edgecolor   {b}                       Line color
%       edgecolor   {[.7 .7 .7]}              Cell line color
%       marker      {.}                       Node marker style
%       markercolor {b}                       Node marker color
%       cellabels   on|{off}                  Show cell numbers
%       nodelabels  on|{off}                  Show vertex/node numbers
%       fontsize    {axes default}            Font size used in text
%       selcells    {all}                     Index vector to cells to plot
%       selsubd     {all}                     Index vector to subdomains to plot
%       renderer    {empty string}|plotly     Render in MATLAB or Plotly
%       plotlyfile  {default}                 Filename to use for plotly export
%       parent      {gca}                     Plot axes handle
%       setaxes     {on}|off                  Use axes modification settings
%
%   See also PLOTBDR, PLOTSUBD, POSTPLOT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help plotgrid, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'plotgrid', varargin{:} );
if( ~nargout ), clear varargout; end
