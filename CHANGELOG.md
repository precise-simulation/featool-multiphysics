 FEATool Multiphysics Changelog
================================


2025-02-14 version 1.17.3
-------------------------

- Warning for MATLAB 2025a GUI incompatibility and performance
- Add checkbox toggle for flow mode automatic pressure contstraint
- Warning for incompatible (openfoam.org foundation) OpenFOAM installation
- Minor bug fixes

2024-11-27 version 1.17.2
-------------------------

- OpenFOAM/SU2 add option to save casefiles on error
- GUI consistency improvements
- SU2 default to V2003m SST turbulence model
- Minor bug fixes


2024-10-20 version 1.17.1
-------------------------

- Added categories to physics mode selection list
- New import fea struct from MATLAB menu option
- Recommendation to use ESI openfoam.com (not openfoam.org)
- Minor bug fixes


2024-09-30 version 1.17
-----------------------

- Physics mode for compressible high Ma number turbulent flows
  + Add tutorial model for supersonic compressible flow past a prism
- Updated OpenFOAM solver interface
  + Exposed OpenFOAM MATLAB API
  + Support import of OpenFOAM dictionaries/data
  + Support coupled flow + heat transfer
  + Support for multiple subdomains
  + Support for chtMultiRegionFoam, buoyantBoussinesqFoam, and sonicFoam solvers
  + CFD script model support for the OpenFOAM solver
- SU2 support for turbulent compressible and high Ma flows
- Add new heat exchanger and heat transfer script model examples
- Improved UI performance and responsiveness


2024-06-25 version 1.16.6
-------------------------

- Preliminary OpenFOAM support for conjugate heat transfer
  with multiple subdomains (chtMultiRegionFoam/SimpleFoam)
- MacOS App OpenFOAM support (https://github.com/gerlero/openfoam-app)
- Script example for conjugate heat transfer (ex_heattransfer10)
- Removed averaging of initial conditions for boundaries shared between subdomains
- Improved performance of Robin boundary condition assembly
- Fix for 3D geometry plotting (identification and deduplication of shared boundaries)


2024-03-15 version 1.16.5
-------------------------

- Support CSG formulas for robust (STL/OBJ) meshes
- New supersonic flow passed a wedge tutorial
- Change supersonic flow over bump tutorial to 3D
 - Support HOME/Documents/MATLAB/.featool/featool.ini configuration
  + performance mex and mumps solver disabled by default
- Fix bug in geometry object rotation
- Minor fixes


2023-11-30 version 1.16.4
-------------------------

- Fix Excel (xlsx) data export
- Add Data Export dialog box
  + Support custom data export expressions
  + Support custom data export coordinates
- Improved app startup and shutdown time


2023-09-10 version 1.16.3
-------------------------

- Add rmphys function and equation mode menu option
  to remove physics modes
- Improved boundary plots for plotsubd
- Updated OpenCASCADE geometry kernel to v7.7.2
- Fix 2D STEP/IGES geometry import
- Fix plot 3D block/hexahedal grid
- Fix for Matlab 2023b


2023-03-12 version 1.16.2
-------------------------

- Show axes coordinate system orientation for 3D views
- Update MUMPS linear solver
  + Support all MATLAB versions (Windows and Linux)
  + Support for Intel MKL BLAS, OpenMP, and >2GB large arrays
- Updated SU2 Code CFD solver to v7.5.0
- Add functionality to show geometry object information


2022-10-20 version 1.16.1
-------------------------

- Updated m-file export output to support the MATLAB "publish" command
- New "Create Model Report" menu option to Generate model reports in
  Html, PDF, Microsoft Word/PDF, Latex, and XML formats
- Updated OpenCASCADE geometry kernel to v7.6.3
- Updated SU2 Code CFD solver to v7.4.0
- Add new surface currents/skin effects tutorial model


2022-08-29 version 1.16
-----------------------

- Performance improvements for built-in multiphysics solver
- Improved 3D geometry rendering performance
- New 3D geometry defeaturing functionality
- Support for edge/vertex selection for chamfer/fillet operations
- Support for PDF image and Excel data export
- Improved save/load model file performance


2022-06-06 version 1.15.5
-------------------------

- Fix for STEP geometry import with >= 2 objects
- Improved and faster expression evaluation in assembly


2022-05-09 version 1.15.4
-------------------------

- Added 2D geometry tool for Bezier and spline curves
- Updated OpenCASCADE geometry kernel to v7.6.2
- Support for binary brep (.bin) geometry format
- Updated SU2 Code CFD solver to v7.3.1
- Added k-Omega wall function support for SU2 solver
- Stand-alone versions for Windows and Linux


2021-12-01 version 1.15.1
-------------------------
- Geometry import option from bitmap image (bmp, jpeg, png)
- Support for structured meshing of geometry primitives
- Updated plotly library to version 2.6.2
- Updated ParaView Glance library to version 4.17.1,
  and added support for slice and isosurface plot options
- Linux support for HiDPI screens
- Mouse controls for 3D zoom, pan, and rotate


2021-11-01 version 1.15
-----------------------

- Added non-Newtonian flow physics mode for power law
  Bird-Carreau, and Cross power law viscosity models
- New non-Newtonian tutorial example for polymer flow extrusion
- New 3D solid-thermal stress validation and tutorial example
- Improved graphics performance for 3D plots
- Changed 3D plots to fill the whole viewport with
  CAD style orbit, pan, and zoom controls
- Added flip/reverse view option by double clicking
  on the 3D xy, xz, and yz quickview buttons
- Added solve mode GUI menu option for using solver hooks
- Improved FEniCS data and file exchange performance


2021-08-30 version 1.14.3
-------------------------

- Various performance improvements
- Added view boundaries/subdomains menu option
  (for specifying mesh sizes on individual geometric entities)
- Added functionality for anisotropic meshing


2021-07-26 version 1.14.2
-------------------------

- Added electro-osmotic flow multiphysics tutorial example
- Fix for STEP geometry export
- Improvments to STL geometry import
- FEniCS support for 3D models with >= 2 subdomains
- OpenFOAM changed to use of consistent symmetryPlane boundary conditions


2021-05-24 version 1.14.1
-------------------------

- Preliminary support for built-in and robust 3D mesh generators
- GUI menu option to manually renumber/reassign 3D boundaries
- Heated pipe with cooling fins validation tutorial


2021-04-19 version 1.14
-----------------------

- Release license menu option to move license to a new computer system
- Raw STL import and Gmsh meshing option for difficult STL geometries
- Solver restart options using partial solution vectors
- Coupled OpenFOAM and FEATool Multiphysics multi-solver simulation model
- Improved and optimized treatment of slip boundary conditions


2021-03-13 version 1.13.3
-------------------------

- Updated OpenFOAM and SU2 interfaces to support parallel computations
- Added OpenFOAM binary GUI option to support user defined FOAM solvers
- Added support for ESI-OpenCFD native Windows OpenFOAM distribution


2021-03-10 version 1.13.2
-------------------------

- Updated FEniCS interface with full support for parallel computations
- Change to binary hdf5 format for FEniCS mesh and solution data
- New impexp_hdf5 and fenics_import utility functions
- Added point/edge/integral source and constraint buttons
  to the subdomain/boundary modes toolbars


2021-02-17 version 1.13.1
-------------------------

- Fixes for mesh generation of point objects
- Updated OpenCASCADE geometry kernel to v7.5.1
- Updated SU2 Code CFD solver to v7.1.0
- Select UMFPACK instead of MUMPS as default linear solver
  for MATLAB 2020a and later due to Mex file incompatibilities


2020-10-20 version 1.13
-----------------------

- New 3D CAD geometry modeling functionality
  - Chamfer (bevel edges)
  - Fillet (round edges)
  - Extrude/Revolve faces to solids
  - Loft faces
  - Split objects by cutline/cutplane
  - Sketching in 2D Workplanes
- Support for user defined solver monitor functions
- View, edit, and modify OpenFOAM case and dictionary files


2020-07-01 version 1.12.4
-------------------------

- Improved GEOMTool geometry engine with planar 2D support,
  faster and more robust geometry decomposition/analysis
  - Added chamfer and fillet geometry operations
  - Added split geometry object by cutplane/cutline
- Algrebraic Multigrid (AMG) linear solver option


2020-05-29 version 1.12.3
-------------------------

- Support for boundary and subdomain integral constraints
- Generalized slice plot postprocessing option
- Gradient reconstruction for evaluation of first order derivatives
- Animation of multiple time steps/solutions postprocessing option
- Streamline plot option
- Updated postprocessing dialog box GUI
- Improved OpenFOAM and SU2 discretization scheme options


2020-04-22 version 1.12.2
-------------------------

- Support for multiple imported CAD geometries
- Improved Gmsh import
- OpenFOAM and SU2 support for pressure inlet boundary conditions
- Added arrow plot scaling option
- Added solver settings GUI option for numerical integration order


2020-04-09 version 1.12.1
-------------------------

- Full support for FEniCS external solver
  - Supports conversion of most FEATool models
    to FEniCS Python simulation scripts
  - Time dependent and nonstationary problems
  - Nonlinear and discontinuous coefficients
- Improved postprocessing and visualization
  - Faster and more efficient 3D graphics
  - GUI support for more colors and colormaps,
    titles, and plot configurations
- GUI functionality to view partial grids
- Faster GUI command/terminal output speed
- Improved and more accurate Robin boundary condition treatment
- Faster expression evaluation via nodal interpolation
- New electrostatics micro strip capacitance tutorial


2020-02-23 version 1.12
-----------------------

- New SU2 CFD solver with built-in GUI and CLI interfaces
- Built-in algebraic turbulence model for Navier-Stokes physics mode
- Support for several OpenFOAM turbulence model variants
  (Spalart-Allmaras, k-epsilon (RNG and realizable), and k-omega (SST))
- New BREP based 3D CAD engine with support for BREP, IGES, OBJ,
  STEP, and STL CAD format import, export, and conversion
- Support for 2D boundary layers in Gmsh


2019-10-14 version 1.11.1
-------------------------

- Support for IGES and STEP CAD file import
- Improved feature and boundary recognition for STL CAD file import
- Geometry and grid export to web browser
- Support for parsing space derivatives in time equation coefficients


2019-08-01 version 1.11
-----------------------

- Dedicated and improved treatment of internal/interior boundaries
- Fluid-structure interaction physics mode and solver
- MUMPS and iterative (GMRES and BiCGStab with ILU(k)) linear solvers
- Support for Gmsh v3/4 GEO file format geometry import and export
- New polyhedron, ellipsoid, cone, and torus geometry object primitives
- Support for external BRL-CAD and OpenCASCADE (Gmsh) geometry engines
- New examples and tutuorials for thermo-mechanical couplings,
  fluid-structure interaction, and non-Newtonian flow
- Added P1/P2 bubble, Nedelec curl, and Raviart-Thomas vector
  finite element shape functions


2019-04-15 version 1.10
-----------------------

- Finite element GUI script (fes) file format and built-in tutorials
- Support for edge constraints in 3D
- Support multiple STL sections and files during STL CAD file import
- Eigenvalue solver
- Functionality to make geometry object copies
- Deformation postprocessing plot option
- Support for spaces in paths/filenames in external grid generator
- Improved evaluation of expressions in grid points
- GUI improvements


2018-11-30 version 1.9
----------------------

- New compressible Euler equations physics mode
  (inviscid compressible flows)
- New swirl flow physics mode (axisymmetric non-zero
  azimuthal velocity)
- OpenFOAM CFD solver integration improvements
  - Support for general expressions (non-constant)
    in initial and boundary conditions
  - Support axisymmetry and swirl flows
  - Support for inviscid compressible flows
  - Improved k-epsilon/omega turbulence model interface
  - Realtime convergence curves plots
- NACA 4-series pre-defined wing geometry object
- Potential flow velocity field initialization
- Improved resolution of curved geometry boundaries
- Improved parametrization and meshing of curved boundaries
- Automatic shock capturing and stabilization for
  convection dominated flow regimes
- Support for first order P1P1 discretization and
  stabilization of flow problems
- Extended backwards compatibility to MATLAB 2009b
- Click and point evaluation for 2D surface plots
- Added built-in GUI tutorials and model examples


2018-05-17 version 1.8
----------------------

- GUI and CLI support for the external OpenFOAM CFD solver
- New CSG BSP tree based geometry engine
  The new geometry engine is used by default but original
  engine may be enabled with a set_geomlib(1) function call
- STL CAD geometry import and export support in 3D and 2D
- Cleaned up and consistent grid generation interfaces
  (identical boundary and subdomain numbering)
- Cached equation and expression parsing for faster evaluations
- Revised nonlinear Newton assembly and solver framework
  with automatic numeric Jacobian assembly calculation


2018-01-08 version 1.7.1
------------------------

- Built-in support for the external mesh generators Gmsh and Triangle
- Improved Gmsh grid import/export functionality with support for
  boundary information
- Redesigned and simplified grid mode toolbar
- GUI stability improvements and fixes
- GUI geometry export in Gmsh geo and Triangle poly formats
- gridextrude: support for extrusion of unstructured triangular
  simplex grids
- Swirl flow and periodic model examples
- Improved sparse matrix modification performance for Dirichlet
  boundary conditions and non-participating dependent variables
- Social sharing of postprocessing images and results


2017-08-10 version 1.7
----------------------

- FEniCS solve mode integration
- Plotly post mode integration
- Import/export of Dolfin XML grid format
- New physics modes:
  - Darcy's Law
  - Brinkman Equations
  - Electrostatics
  - Euler-Bernoulli Beam
  - Magnetostatics
- Built-in support for axisymmetry
- Solver support for active/inactive dependent variables
  in different subdomains and boundary couplings
- Added function to reconstruct inner boundaries
- Explicit construction of geometry object
  distance functions
- Support for 2nd order derivatives
  and equation definitions
- Parsing of groups of derivatives in parentheses
- New models and examples
- Support arbitrary order quadrature rules for
  all cell types and renamed cubrule* to quadrule*


2017-01-23 version 1.6
----------------------

- Support for external CFD solvers
- 3D fluid flow model examples using external CFD solvers
- Automatic quadrilateral grid generation in 2D
- Uniform grid refinement in 3D
- Projection to geometrical boundaries for uniform
  2D grid refinement and grid cell conversion
- 3 new axisymmetric stress-strain model examples
- New non-linear PDE models (in 1D and 2D ),
  and one new 2D Laplace equation test model
- 3rd-5th order Lagrange shape functions
- Cubic Hermite shape functions in 1D and 2D
- Support for tetrahedral quadrature up to order 5


2016-06-07 version 1.5
----------------------

- Redesigned 3D geometry engine
  - Calculates geometry boundary faces and edges for
    more precision in 3D geometry and grid generation
- 3D slice plot refactored and improved for new geometry engine
- New conforming hex to tet grid conversion method (1 Hex -> 28 Tets)
- GUI menu support for grid conversion and smoothing
- Allow partial rings in ringgrid generation function
- Generalized Navier-Stokes slip boundary conditions
  (supports curved and non-axis aligned boundaries)
- Allow external function calls for coefficients
- Progress waitbars for grid generation and solvers


2016-01-28 version 1.4
----------------------

- New physics mode features:
  - Heat transfer natural convection and radiation boundary conditions
  - Temperature stress-strain coupling in structural mechanics physics modes
  - Artificial and (anisotropic) streamline diffusion for stabilization
    of convection dominated high Re/Pe flows
  - Slip boundary condition option for the Navier-Stokes equations
- Grid import and export in GiD, Gmsh, and Triangle formats
- Define point sources and constraints
- New postprocessing features:
  - Save and export plots to image files (jpeg, png, eps)
  - Surface height plot option in two dimensions
  - Max and min boundary and subdomain evaluation functions
  - GUI option to evaluate general expressions on points and lines
  - Set colorbar min and max limits
  - Postprocessing plots show a title describing the plot
- New models and examples
- Solver, core, and grid generation fixes and enhancements


2015-08-30 version 1.3
----------------------
- Improved 2D geometry engine (especially boundary assignment
  and treatment of multiple geometry objects)
- Fractional step-theta scheme and HRZ mass lumping for
  time dependent problems
- Newton option for nonlinear stationary problems
- Optimized unstructured grid generation with better treatment
  of multiple subdomains
- Saving to m-script file
- Postprocessing and visualization export to Plotly (web)
- GUI option for model constants and expressions


2015-03-03 version 1.2
----------------------
- Full 3D support
- Implicit (distance/levelset function based) geometry modeling engine
- Unstructured grid generation
- Strong to weak equation formulation parser
- Custom equation(s) physics mode
- Extensibility, call to external grid, solver, and postprocessing supported
- Optimized FEM assembly and transformation Jacobian calculation routines
- Minor bug fixes


2014-11-27 version 1.1
----------------------
- Minor bug fixes


2014-05-17 version 1.0
----------------------
- Initial release


2013-07-25 version 0.9
----------------------
- Preliminary version
