%IMPEXP_GMSH Import/export grid in GMSH format.
%
%   [ GRID ] = IMPEXP_GMSH( FILE_NAME, MODE, DATA, IIND, FID_LOG, BIX )
%   Import and export of Gmsh (.msh) Ascii mesh and grid file format.
%
%   FILE_NAME is a string specifying the (root) file name to process.
%
%   MODE can either be a string indicating IMPORT (where boundary
%   numbering are reconstructed from Gmsh physical or elemental tags if
%   available), IMPORT_BDR (boundary reconstruction with gridbdr), or
%   EXPORT.
%
%   With the (default) IMPORT mode option the function returns a GRID
%   struct. While for export, the input DATA can be given as either
%   a valid fea finite element data struct, or a FEATool grid struct.
%
%   IIND is an optional array for physical domains (subdomains in 2D and
%   boundaries in 3D) to extract and import (default all).
%
%   FID_LOG is an optional log file handle for message output
%   (negative for Gui output or empty for no output).
%
%   BIX is an optional index vector to correct boundary
%   numbering. Alternatively, BIX = 0/false disables boundary
%   reconstruction during import.

% Copyright 2013-2024 Precise Simulation, Ltd.
