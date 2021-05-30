function [ varargout ] = plotsubd( varargin )
%PLOTSUBD Subdomain plotting function.
%
%   [ H ] = PLOTSUBD( PROB, VARARGIN ) Function to plot and highlight
%   subdomains. PROB is a valid finite element problem struct.
%   Accepts the following property/value pairs.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       selsubd     {all}                     Index vector to subdomains to plot
%       axis        off|{on}                  Show axes
%       grid        on |{off}                 Show grid
%       normals     on |{off}                 Show normals
%       axequal     off|{on}                  Axis equal setting
%       bbox        {0.05}                    Size of bounding box (0=off)
%       labels      on/{off}                  Print boundary numbers
%       fontsize    {2*axes default}          Font size used in text labels
%       view        {[-35 20]}                3D view setting
%       alpha       {0.3}                     Transparency level
%       parent      {gca}                     Axes handle to plot in
%
%   See also PLOTGRID, PLOTBDR, POSTPLOT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help plotsubd, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'plotsubd', varargin{:} );
if( ~nargout ), clear varargout; end
