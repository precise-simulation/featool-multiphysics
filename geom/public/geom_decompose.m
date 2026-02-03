%GEOM_DECOMPOSE Analyze and decompose geometry.
%
%   [P,E,F,V,IP] = GEOM_DECOMPOSE( GEOM, TOL, IS_ANALYZED ) Decomposes
%   the geometry struct GEOM into unique points/vertices P, edges/line
%   segments E, faces or surfaces F, and volumes V (for 3D geometries).
%   IS_ANALYZED (default false) indicates if the geometry has already
%   been processed by GEOM_ANALYZE. The optional input argument TOL
%   defines the tolerance for deduplication and zeroing vertex
%   coordinates (default eps*1e3). IP indicates which vertices in P
%   are due to point geometry objects (and not part of edges, faces,
%   or volumes).
%
%   Geometries consisting of multiple geometry objects will be analyzed
%   and if required split into separate non-overlapping objects. The
%   decomposition process tries to apply one intersection and two
%   subtraction operations to each combination of two geometry separate
%   objects.
%
%   The output argument P is a array of unique points/vertices of
%   size (n_p,n_sdim).
%
%   The two last columns of the edge/line segment array E (size (n_e,3))
%   contains start and end indices to vertices in P that make up the
%   edges. The first column indicates edge segment groupings
%   (corresponding to boundaries in 2D and geometry edges in 3D,
%   No defined edge grouping is indicated with zeros in 3D).
%
%   The face/surface array F of size (n_f,2+n_max_e) consists of
%   indices to edges that make up the faces (a negative entry indicates
%   reversed edge direction). Similar to the first column of E, the first
%   column of F also indicates face groupings (integer > 0). The second
%   boolean column indicates if each face is an outer boundary (0) or a
%   hole (1). For faces with less than n_max_edges (the maximum number
%   edges in any face) the extra entries are padded with zeros.
%
%   The volume array V (for 3D) indicates which faces F constitutes the
%   volumes (size n_v,2+n_max_f) where a negative entry also here
%   indicates reversed face direction. Similarly to F, the first column
%   specifies the volume group, and the second if the entry volume
%   describes the outer (external) or inner volume boundary (hole).
%
%   See also GOBJ_DECOMPOSE, GEOM_ANALYZE

% Copyright 2013-2026 Precise Simulation, Ltd.
