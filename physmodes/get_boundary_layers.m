function [ varargout ] = get_boundary_layers( varargin )
%GET_BOUNDARY_LAYERS Find boundary layers.
%
%   [ BL ] = GET_BOUNDARY_LAYERS( FEA, IND_BL, HL, NL, AR, ARG )
%   Find boundary layers in FEA struct and return the struct BL
%   with the following fields for connected wall boundaries
%
%       boundaries - list of boundaries in layer
%       svertices  - list of single vertices in layer
%       sedges     - list of single edges in layer (only 3D)
%       gmsh       - gmsh boundary field
%
%   The optional array IND_BL is used to override automatic detection
%   of wall boundaries (boundary condition selection 1 for
%   Navier-Stokes, swirl flow, Brinkman equations, and fluid-structure
%   interation physics modes). If IND_BL is empty layers for all
%   boundaries will be returned.
%
%   The optional numerical arguments HL, NL (default 5), and AR
%   (default 1.2) specifies the layer height, number of cells in
%   layer, and aspect ratios for gmsh output.
%
%   ARG is an optional struct to give parameters for Gmsh boundary
%   layer generation with the fields hmax, hmaxb, and thmax (default 5*pi/3).

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help get_boundary_layers, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'get_boundary_layers', varargin{:} );
if( ~nargout ), clear varargout; end
