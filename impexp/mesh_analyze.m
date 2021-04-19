function [ varargout ] = mesh_analyze( varargin )
%MESH_ANALYZE 3D surface mesh analysis and feature recognition.
%
%   [ IND_FG, IND_EG, H_PLOT ] = MESH_ANALYZE( P, F, N, TOL_TH,
%   TOL_AR, TOL_CO, OPT, I_PLOT, H_AX ) Perform analysis and
%   feature recognition on the STL surface mesh data given by the
%   vertex coordinates in P (n_p x 3), faces F (n_f x 3), and outward
%   directed normal vectors N (n_f x 3).
%
%   TOL_TH is an optional vector with four entries specifying the
%   following angles (in degrees); dihedral angle cutoff to classify
%   sharp features (default 30), angles to classify needle and cap
%   facets for identifying non-sharp features (default 10 and 160),
%   and maximum edge deviation for growing unconnected feature edges
%   (default 20).
%
%   TOL_AR optionally specifies various facet related aspect ratio
%   related tolerances for non-sharp feature identification, namely;
%   Minimum facet height to length (hl-)ratios (default 3), max to min
%   height-ratio (default 2.5) for keeping non-sharp features related to
%   cap edges, as well as hl-ratio for eliminating single unconnected
%   edges (default 1.2).
%
%   Coalescence/merging of small surface patches with less than
%   TOL_CO(1) facets (default 2) and maximum edge TOL_CO(2) angle
%   deviation (default 30 degrees) is also optionally performed.
%
%   The additional flags in OPT are SHARP and EXTEND. SHARP toggles
%   splitting of sharp (>100 degrees) edges which helps avoid complex
%   shaped patches. SHARP=2 (default) adds all potential edges around
%   the sharp angled vertex, and SHARP=1 only adds the edge with least
%   angle deviation. The flag EXTEND toggles extending/growing of
%   disconnected edges (default true).
%
%   The parameter I_PLOT enables plot generation, 1 for
%   surface/boundary plot, and 2 for edge feature plot. An axis handle
%   H_AX can optionally be give.
%
%   Output arguments are the facet and edge groupings in IND_FG and
%   IND_EG (corresponding to edges computed with GRIDBDRE, and where
%   IND_EG is zero for internal facet edges). And H_PLOT contains any
%   graphics handles.
%
%   See also MESH_REPAIR, GRIDBDRE

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help mesh_analyze, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'mesh_analyze', varargin{:} );
if( ~nargout ), clear varargout; end
