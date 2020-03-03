function [ varargout ] = gridgen( varargin )
%GRIDGEN Grid generation for geometry objects.
%
%   [ GRID, STATS ] = GRIDGEN( SIN, VARARGIN ) Function to generate a
%   grid for geometry objects by calling an external grid generation
%   algorithm. SIN is a valid fea problem struct or cell array or
%   geometry objects. Accepts the following property/value pairs.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       gridgen     string {default}          Grid generation algorithm: default,
%                                                gmsh, gridgen2d, mesh2d, triangle, or quad
%       hmax        scalar/arr {0.1}          Target grid size for subdomains
%       hmaxb       scalar/arr {[]}           Target grid size for boundaries
%       avhb        logical {true}            Average hmax to boundaries
%       fixpnt      array   [n_p,n_sdim]      Array of (fixed) points that must be present
%       nsm         scalar  {3}               Number of (post) grid smoothing steps
%       intb        logical {true}            Output interior/internal boundaries
%       waitbar     scalar  {0}               Show/hide waitbar
%       fid         scalar  {1}               File identifier for output ([]=no output)
%
%   GRIDGEN specifies which grid generation algorithm to use and calls
%   the corresponding grid generation code.
%
%   HMAX indicates target grid cell diameters, and is either a numeric
%   scalar prescribing the grid size for the entire geometry, or an
%   array with HMAX values corresponding to individual
%   subdomains. Positive HMAX values uses the minimum mesh size for
%   shared boundaries, while negative applys the mean value.
%
%   HMAXB is analogous to HMAX but related to boundaries
%   (edges/faces). HMAXB can be a single scalar applicable to all
%   boundaries, a numeric array with entries corresponding to
%   individual boundaries.
%
%   The AVHB logical flag toggles if internal boundaries inheriting
%   HMAX values (when HMAXB is unspecified) should be assigned the
%   smallest HMAX value from neighbouring subdomains, or the mean value
%   (default).
%
%   NSM (default 3) the number of GRIDSMOOTH smoothing steps to perform.
%
%   INTB toggels interior/internal boundaries on (default) or off.
%
%   Examples:
%
%      1) Unit square with uniform global grid size set to 0.1.
%
%      grid = gridgen( {gobj_rectangle()}, 'hmax', 0.1 );
%      plotgrid( grid )
%
%      2) Unit square with a finer grid along the top boundary (using Triangle).
%
%      grid = gridgen( {gobj_rectangle()}, 'hmax', 0.5, ...
%                      'hmaxb', [0.5 0.5 0.01 0.5], 'gridgen', 'triangle' );
%      plotgrid( grid )
%
%      3) Domain with curved boundaries meshed with quadrilaterals (using Gmsh).
%
%      geom.objects = {gobj_rectangle() gobj_circle([0 0],.6) gobj_circle([1 1],.3,'C2')};
%      geom = geom_apply_formula( geom, 'R1-C1-C2' );
%      grid = gridgen( geom, 'hmax', 0.1, 'gridgen', 'gmsh', 'quad', true );
%      plotgrid( grid )
%
%      4) Two connected subdomains with a shared boundary (using Gridgen2D).
%
%      geom.objects = { gobj_polygon([-2e-3 -8e-3;0 -8e-3;0 -6e-3;0 6e-3;0 8e-3;-2e-3 8e-3]), ...
%      gobj_polygon([0 -6e-3;2e-3 -5e-3;2e-3 4e-3;0 6e-3]) };
%      hmax  = 5e-4;
%      hmaxb = hmax*ones(1,4);
%      hmaxb(9) = hmax/5;
%      grid  = gridgen( geom, 'hmax', hmax, 'hmaxb', hmaxb, 'gridgen', 'gridgen2d' );
%      plotgrid( grid )
%
%      5) Composite component with two subdomains (using Gridgen2D).
%
%      r1 = gobj_rectangle( 0, 0.11, 0, 0.12,  'R1' );
%      c1 = gobj_circle( [ 0.065 0 ],   0.015, 'C1' );
%      c2 = gobj_circle( [ 0.11 0.12 ], 0.035, 'C2' );
%      c3 = gobj_circle( [ 0 0.06 ],    0.025, 'C3' );
%      r2 = gobj_rectangle( 0.065, 0.16, 0.05, 0.07, 'R2' );
%      c4 = gobj_circle( [ 0.065 0.06 ], 0.01, 'C4' );
%      geom.objects = { r1 c1 c2 c3 r2 c4 };
%      geom = geom_apply_formula( geom, 'R1-C1-C2-C3' );
%      geom = geom_apply_formula( geom, 'R2+C4' );
%
%      grid  = gridgen( geom, 'hmax', [0.0025 0.01 0.0025], 'gridgen', 'gridgen2d' );
%      plotgrid( grid )
%
%      6) Complex geometry with several holes and subdomains (using Gridgen2D).
%
%      w = 10e-4; L = 3*w; H = 5*w;
%      p1  = gobj_polygon( [w/10 0;(L-w/4)/2 0;(L-w/4)/2 H;0 H;0 H/3], 'P1' );
%      p2  = gobj_polygon( [(L+w/4)/2 0;L 0;L H-H/3;L H;(L+w/4)/2 H], 'P2' );
%      r1  = gobj_rectangle( (L-w/4)/2, (L+w/4)/2, 0, H, 'R1' );
%      c1  = gobj_circle( [2*w/3 3*w], w/3, 'C1' );
%      c2  = gobj_circle( [2*w/3 2*w], w/3, 'C2' );
%      c3  = gobj_circle( [2*w/3 1*w], w/3, 'C3' );
%      c4  = gobj_circle( [L-w/2 4.5*w], w/8, 'C4' );
%      c5  = gobj_circle( [L-w   4.5*w], w/8, 'C5' );
%      c6  = gobj_circle( [L-w/2 4*w], w/8, 'C6' );
%      c7  = gobj_circle( [L-w   4*w], w/8, 'C7' );
%      c8  = gobj_circle( [L-w/2 3.5*w], w/8, 'C8' );
%      c9  = gobj_circle( [L-w   3.5*w], w/8, 'C9' );
%      c10 = gobj_circle( [L-w/2 3*w], w/8, 'C10' );
%      c11 = gobj_circle( [L-w   3*w], w/8, 'C11' );
%
%      geom.objects = { p1 p2 r1 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 };
%      geom = geom_apply_formula( geom, 'P1-C1-C2-C3' );
%      geom = geom_apply_formula( geom, 'P2-C4-C5-C6-C7-C8-C9-C10-C11' );
%
%      hmaxb = zeros(1,21);
%      hmaxb([6 20]) = w/50;
%      grid = gridgen( geom, 'hmax', w./[5 5 20], 'hmaxb', hmaxb, 'gridgen', 'gridgen2d' );
%      plotgrid( grid )

% Copyright 2013-2020 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridgen, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridgen', varargin{:} );
if( ~nargout ), clear varargout; end
