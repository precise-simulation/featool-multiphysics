function [ varargout ] = impexp_stl( varargin )
%IMPEXP_STL Import/export geometry or grid in STL format.
%
%   [ OUT ] = IMPEXP_STL( FILE_NAME, MODE, DATA, VARARGIN ) Import
%   and export of CAD geometries and meshes in STL format. Both 3D and
%   planar 2D geometries are supported (planar STL data will be
%   converted to 2D geometry objects).
%
%   FILE_NAME can either be a string to specify the file name to
%   process, or a cell array of strings specifying several STL
%   files. where in each file can be parsed in to separate subdomains
%   or boundary regions.
%
%   MODE is a string indicating either IMPORT or EXPORT. For the
%   default IMPORT mode option the function returns a FEATool geometry
%   struct, grid struct, or unprocessed STL data. For export a postfix
%   string (for example export_R1) can be used to assign a specific
%   geometry object to export (and not split 3D boundaries to solid
%   sections). The RAW flag can be used to return an unprocessed
%   geometry object, or the STL vertices and facets.
%
%   During export, DATA can either be a FEATool finite element
%   geometry or mesh struct to export to STL format. For 3D import
%   DATA is an optional geometry struct where the objects are used to
%   assist with assigning and matching facets to boundaries (instead
%   of using solid sections or mesh analysis to group facets). An
%   additional optional logical field nrev in the geometry object
%   structs is supported for reversing of the reference normals (due
%   to subtraction csg operation).
%
%   Accepts the following optional property/value pairs during import.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       solid       scalar  {1}               Solid section processing
%                                                0 - Merge solid sections and analyze mesh
%                                                1 - Perform mesh analysis without merging
%                                                2 - Assign solid sections as boundaries
%       fsubd       logical {false}           Process files as separate subdomains
%                                             Mesh analysis parameters (see MESH_ANALYZE)
%       tol_th      array   {[30,10,160,20]}     Angle tolerances for feature detection
%       tol_ar      array   {[3,1,1.5]}          Aspect ratio tolerances
%       tol_co      array   {[2,30]}             Path coalescence parameters
%       sharp       integer {2}               Split regions with sharp edges
%       extend      logical {true}            Extend/grow unconnected edges
%       raw         integer {0}               1 = Import raw STL geometry object,
%                                             2 = Import raw STL data (vertices, facets)
%       check       logical {true}            Check boundary connectivities
%
%   See also MESH_REPAIR, MESH_ANALYZE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help impexp_stl, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'impexp_stl', varargin{:} );
if( ~nargout ), clear varargout; end
