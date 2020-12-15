function [ varargout ] = plotedg( varargin )
%PLOTEDG Boundary edges plot function.
%
%   [ H ] = PLOTEDG( PROB, VARARGIN ) Plot and highlight 3D boundary
%   edges. PROB is a valid finite element problem structure. Accepts
%   the following property/value pairs.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       axis        off|{on}                  Show axes
%       grid        on |{off}                 Show grid
%       bbox        {0.05}                    Size of bounding box (0=off)
%       linewidth   {2}                       Line width
%       colors      {r,g,b,c,m,y}             Cell array with colors
%       labels      on/{off}                  Print boundary numbers
%       fontsize    {2*axes default}          Font size used in text labels
%       view        {[-35 20]}                3D view setting
%       alpha       {0.5}                     Transparency level
%       seledg      {all}                     Index vector to edges to plot
%       parent      {gca}                     Axes handle to plot in
%
%   See also PLOTGEOM, PLOTGRID, PLOTBDR, PLOTSUBD, POSTPLOT

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help plotedg, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'plotedg', varargin{:} );
if( ~nargout ), clear varargout; end
