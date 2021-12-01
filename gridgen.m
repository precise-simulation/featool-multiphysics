function [ varargout ] = gridgen( varargin )
%GRIDGEN Grid and mesh generation for geometry objects.
%
%   [ GRID, STATS ] = GRIDGEN( SIN, VARARGIN ) Generates a grid/mesh
%   for the geometry defined by the objects in SIN by calling an
%   external grid generation algorithm. SIN can be either a valid fea
%   problem struct with geometry defined in SIN.geom.objects, a cell
%   array of multiple geometry objects, or a single geometry
%   object. Accepts the following property/value pairs are available
%   with the default grid generation algorithm.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       gridgen     string {default}          Grid generation algorithm: default, gmsh
%                                                robust (3D), gridgen2d, or triangle (2D)
%       dprim       logical    {true}         Structrured meshing of primitives
%       hmax        scalar/arr {0.1}          Target grid size for subdomains
%       hmaxb       scalar/arr {[]}           Target grid size for boundaries
%       hmaxe       scalar/arr {[]}           Target grid size for edges (3D)
%       grading     scalar     0.3            Mesh grading/growth rate
%       eledge      scalar     1              Elements per edge
%       quad        logical    {false}        Use quad meshing (2D)
%       intb        logical    {true}         Output interior/internal boundaries
%       waitbar     scalar     {0}            Show/hide waitbar
%       fid         scalar     {1}            File identifier for output ([]=no output)
%
%   GRIDGEN specifies which grid generation algorithm to use and calls
%   the corresponding grid generation code (default, gmsh, robust
%   (3D), gridgen2d (2D), or triangle (2D)). Enabling DPRIM with the
%   default grid generator enables structured meshing of single
%   geometry object primitives (and no differing boundary/edge mesh
%   sizing).
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
%   HMAXE is analogous to HMAXB but related to 3D edges (only
%   applicable to default grid generator for 3D geometries).
%
%   GRADING specifies the rate at which smaller cells will grow to
%   larger (0 - 1), and ELEDGE optionally prescribes the minimum
%   elements per edge.
%
%   INTB toggels interior/internal boundaries on (default) or off.
%
%
%   The following additional property/value pairs are available with
%   and specific to the GMSH mesh generator.
%
%       Property    Value/{Default}            Description
%       -----------------------------------------------------------------------------------
%       hmaxp       scal/arr {[]}              Target grid size for vertices
%       nref        scalar   {0}               Number of uniform grid refinements
%       algo2       scalar   {2}               2D mesh algorithm (1=MeshAdapt, 2=Automatic,
%                                              5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
%       algo3       scalar   {1}               3D mesh algorithm (1=Del, 2=New Del, 4=Front
%                                              5=Front Del, 6=Front Hex, 7=MMG3D, 9=R-tree)
%       blayer      struct   {[]}              Data struct for Gmsh 2D boundary layers
%       quad        logical  {false}           Use quad meshing (for 2D)
%       intb        logical  {true}            Output interior/internal boundaries
%       avhb        logical  {true}            Average hmax to boundaries
%       nsm         scalar   {3}               Number of (post) grid smoothing steps
%       palign      scalar   {eps}             Distance tolerance to align point objects
%       tol         scalar   {eps*1e3}         Deduplication tolerance
%       compound    logical  {true}            Use Gmsh compound boundaries
%       mshopt      cell     {}                Cell array of additional Gmsh options
%       mshall      logical  {true}            Output/save all meshed entities
%       mshver      integer  {2}               Gmsh msh file version (1/2/4)
%       verbosity   integer  {5}               Gmsh verbosity/output level
%       syscmd      string   {'default'}       Gmsh system call command
%                                              (default 'gmsh fdir/fname.geo -')
%       fname       string   {'featool_gmsh_UID'}  Gmsh imp/exp file name (root)
%       fdir        string   {tempdir}         Directory to write help files
%       clean       logical  {true}            Delete (clean) Gmsh help files
%
%   NREF (default 0) the number of post uniform grid refinement steps.
%
%   ALGO2 and ALGO3 the Gmsh 2D and 3D mesh generation algorithms.
%
%   QUAD (default 0) toggles Blossom-Quad conversion for 2D geometries.
%
%   The BLAYER flag enables boundary layers for boundaries. May
%   optionally given as array of structs where each entry corresponds
%   to a boundary layer field entry.
%
%   The AVHB logical flag toggles if internal boundaries inheriting
%   HMAX values (when HMAXB is unspecified) should be assigned the
%   smallest HMAX value from neighbouring subdomains, or the mean value
%   (default).
%
%   NSM (default 3) the number of GRIDSMOOTH smoothing steps to perform.
%
%   PALIGN sets a minimum distance over which to re-align mesh
%   vertices to point objects.
%
%   Additional Gmsh options can be provided with the cell array
%   MSHOPT. For example to set the "CharacteristicLengthMax" and "AnisoMax"
%   Gmsh options, MSHOPT could be given as
%
%       {{'Mesh', 'CharacteristicLengthMax', '1'}, {'Mesh', 'AnisoMax', '10'}}
%
%
%   The following additional property/value pairs are available with
%   and specific to the GRIDGEN2D mesh generator.
%
%       Property    Value/{Default}            Description
%       -----------------------------------------------------------------------------------
%       q           scalar   {0.65}            Target quality
%       blayer      logical  {false}           Enable boundary layers
%       avhb        logical  {true}            Average hmax to boundaries
%       nsm         scalar   {3}               Number of (post) grid smoothing steps
%       palign      scalar   {eps}             Distance tolerance to align point objects
%
%   Q (default 0.65) specifies a quality target [0 < q < 1.0].
%
%   The BLAYER flag enables boundary layers for internal boundaries
%   (holes). May optionally given as a vector where the 1st entry
%   specifies relative thickness of stretched layer.
%
%   PALIGN sets a minimum distance over which to re-align mesh
%   vertices to point objects.
%
%
%   The following additional property/value pairs are available with
%   and specific to the TRIANGLE mesh generator.
%
%       Property    Value/{Default}            Description
%       -----------------------------------------------------------------------------------
%       q           scalar   {28}              Minimum target angle (quality)
%       syscmd      string   {'default'}       Triangle system call command (default
%                                              'triangle -I -q%f -j -e -a -A %s.poly')
%       fname       string   {'featool_tri_UID'}  Triangle imp/exp file name (root)
%       fdir        string   {tempdir}         Directory to write help files
%       clean       logical  {true}            Delete (clean) Triangle help files
%
%   Q (default 28 degrees) specifies a minimum target angle (values less than 33 are
%   generally acceptable, while higher values might prevent Triangle convergence).
%
%
%   Examples:
%
%      1) Unit circle with uniform global grid size set to 0.1.
%
%      grid = gridgen( gobj_cone, 'hmax', 0.1 );
%      plotgrid( grid )
%
%      2) Unit square with a finer grid along the top boundary (using Triangle).
%
%      grid = gridgen( gobj_rectangle, 'hmax', 0.5, ...
%                      'hmaxb', [0.5 0.5 0.01 0.5], 'gridgen', 'triangle' );
%      plotgrid( grid )
%
%      3) Domain with curved boundaries meshed with quadrilaterals (using Gmsh).
%
%      geom.objects = {gobj_rectangle() gobj_circle([0 0],.6) gobj_circle([1 1],.3,'C2')};
%      geom = geom_apply_formula( geom, 'R1-C1-C2' );
%      grid = gridgen( geom, 'hmax', 0.1, 'gridgen', 'gmsh', 'quad', true, 'verbosity', 3 );
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
%      hmax = w./[5 5 20];     % Set finer mesh size in subdomain 3.
%      hmaxb = zeros(1,21);
%      hmaxb([6 21]) = w/50;   % Set finer mesh size on the in and outlets, boundaries 6 and 21.
%      grid = gridgen( geom, 'hmax', hmax, 'hmaxb', hmaxb, 'gridgen', 'gridgen2d' );
%      plotgrid( grid )
%
%      7) Cone showing the interior y>0 and z<0 (using Gmsh).
%
%      grid = gridgen( gobj_cone, 'hmax', 0.25, 'gridgen', 'gmsh' );
%      plotgrid( grid, 'facecolor', 'none', 'edgecolor', [.7 .7 .7] )
%      plotgrid( grid, 'selcells', '(y>0)&(z<0)' )

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridgen, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridgen', varargin{:} );
if( ~nargout ), clear varargout; end
