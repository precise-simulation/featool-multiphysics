function [ varargout ] = fenics_import( varargin )
%FENICS_IMPORT Import data from FEniCS HDF5 file.
%
%   [ PROB ] = FENICS_IMPORT( FILE_NAME, PROB ) Imports solution and
%   vector data from a FEniCS HDF5 file FILE_NAME, converts it and
%   returns FEATool format. Searches for hdf5 data with FiniteElement
%   attribute and maps it to PROB.SOL.U and PROB.SOL.T if also the
%   timestamp attribute is set in the vector data.
%
%   The data input file must either contain valid hdf5 /mesh and
%   /boundary datasets (as in IMPEXP_HDF5), or optionally a FEATool
%   problem struct PROB with a valid matching grid can be given as an
%   additional input argument.
%
%   Example:
%
%      1) Solve and import a Poisson problem with FEniCS.
%
%      fea.sdim = {'x', 'y'};   % Define and set up model problem.
%      fea.grid = quad2tri(holegrid());
%      fea = addphys( fea, 'poisson' );
%      fea = parsephys(fea);
%      fea = parseprob(fea);
%
%      fea.sol.u = solvestat( fea );   % Solve problem with FEATool Multiphysics.
%
%      fenics( fea, 'modes', {'export','solve'}, ...
%              'fname', 'featool-fenics-test', ...
%              'fdir', pwd(), 'clean', false )   % Export and solve problem with FEniCS.
%
%      fea_1 = fenics_import( 'featool-fenics-test.h5', fea );   % Import FEniCS solution and merge with fea data struct.
%
%      fea_2 = fenics_import( 'featool-fenics-test.h5' );   % Import FEniCS solution and mesh, and reconstruct fea data struct.
%
%      subplot(1,3,1), postplot( fea,   'surfexpr', 'u' )   % Postprocessing and solution comparison.
%      subplot(1,3,2), postplot( fea_1, 'surfexpr', 'u' )
%      subplot(1,3,3), postplot( fea_2, 'surfexpr', 'u' )
%
%   See also FENICS, IMPEXP_HDF5

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help fenics_import, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'fenics_import', varargin{:} );
if( ~nargout ), clear varargout; end
