function [ varargout ] = geom2geo( varargin )
%GEOM2GEO Generate GEO geometry data for use with Gmsh.
%
%   [ GEO ] = GEOM2GEO( SIN, HMAX, HMAXB, HMAXP, FILE_NAME, FID_LOG, TOL,
%                       F_COMPOUND, I_VER, IND_GOBJ, OPS, P_FIX, USE_CAD )
%
%   Converts a geom struct to GEO data for use with the mesh generator Gmsh.
%
%   SIN can either be a valid geometry or precomputed GEO struct, HMAX,
%   and HMAXB are scalars or arrays indicating grid sizes for subdomains
%   and boundaries. FILE_NAME is optional and enables output to file
%   (writes a geo file, default featool_gmsh). FID_LOG specifies a message
%   log file handle (negative for gui output or empty for no output). The
%   optional input argument TOL defines the tolerance for deduplication
%   and zeroing vertex coordinates (default eps*1e3).
%
%   HMAX, HMAXB, and HMAXP are used to calculate and prescribe the
%   mesh sizes. HMAX is a scalar or array prescribing mesh sizes for
%   subdomains, while HMAXB and HMAXP corresponds to boundary and
%   vertex mesh sizes. Mesh sizes are prescribed in Gmsh Restrict
%   Fields for edges and faces, respectively. A non-zero and non-empty
%   HMAX assigns mesh size values to boundary fields connected to
%   subdomains. If HMAX is an array with positive values the minimum
%   value will be used, while negative HMAX values prescribes the mean
%   absolute for each boundary. If HMAXB is given the corresponding
%   boundary values are overwritten (HMAXB overrides HMAX).
%
%   The F_COMPOUND flag indicates if Gmsh compound lines and surfaces
%   should be used (default true). The I_VER interger determines the
%   Gmsh geo output syntax, vertion >= 4 (default) or older v3 syntax.
%
%   IND_GOBJ is an index vector to the geometry objects to output, or
%   empty if the geometry should be analyzed and decomposed in to
%   minimal regions (default).
%
%   The BIX output GEO struct field contains boundary numbers matching
%   to physical followed by additional compound line entries in 2D and
%   surfaces in 3D. BIX can thus be used to set geometrical boundary
%   numbers during import.
%
%   P_FIX is an optional array of points to which should be in the
%   resulting mesh (in addition to point objects).
%
%   USE_CAD (default true) is an optional logical switch which selects
%   use of CAD format files (BREP, STEP, IGES) in Gmsh instead of
%   converting to Gmsh object primitives.
%
%   See also GEO2GEOM

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help geom2geo, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'geom2geo', varargin{:} );
if( ~nargout ), clear varargout; end
