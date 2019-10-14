FEATool Multiphysics™ - _Physics Simulation Made Easy_
======================================================

![FEATool Multiphysics Screenshot](https://raw.githubusercontent.com/precisesimulation/featool-multiphysics/master/featool-multiphysics-screenshot.png)


About
-----

**FEATool Multiphysics™** (short for Finite Element Analysis Toolbox),
is a toolbox for modeling and simulation of coupled physical
phenomena, [partial differential equations](https://en.wikipedia.org/wiki/Partial_differential_equation)
(PDE), continuum mechanics and engineering problems with the
[finite element method](https://en.wikipedia.org/wiki/Finite_element_method)
(FEM).

_FEATool_ aims to provide an easy to use and comprehensive all-in-one
integrated simulation package for all types of finite element and
multiphysics analyses, and combine the best of ease of use, powerful
functionality, and extensibility for both beginners and advanced
users. Features such as an intuitive and easy to learn graphical user
interface (GUI), a complete library of grid generation, and
postprocessing functions, as well as command line interface (CLI)
programming, and interactive and interpreted programming and scripting
capabilities, makes _FEATool_ suitable for everyone from students
learning mathematical modeling, to professionals and engineers wishing
to explore new ideas in a simple, fast, and convenient way.


Key Features
------------

- Easy to use Graphical User Interface (GUI)
- Fully integrated and built-in pre-processing, CAD tools, automatic
  mesh generation, solvers, post-processing, and visualization
  functionality
- Pre-defined physics modes and equations for heat and mass transfer,
  fluid dynamics, structural mechanics, electromagnetics, and
  classical PDE
- Support for custom user-defined equations and PDEs
- Built-in expression parser (enter equations and coefficients _just
  as writing on paper_ without any programming)
- Easy to define fully coupled multiphysics problems
- One-click interfaces for the OpenFOAM®, FEniCS, and Firedrake
  external solvers
- Save and export models in binary file format or editable script
  files (every GUI operation is recorded to a corresponding MATLAB
  function call)
- Fully scriptable and programmable with the MATLAB scripting language
  (including support for integration and embedding in custom
  applications and toolboxes)


System Requirements
-------------------

The _FEATool Multiphysics_ Toolbox can either be used stand-alone as a
fully integrated simulation GUI, or together with MATLAB enabling
command line (CLI) use and scripting. In order to use the _FEATool
Multiphysics_ toolbox, the program files must first be installed on the
intended computer system. Please follow the procedure below for the
corresponding stand-alone or MATLAB toolbox installation modes. Note
that it is recommended to first uninstall older versions of _FEATool_
before updating to a newer version.


Stand-alone Toolbox
-------------------

FEATool Multiphysics has been compiled to run in stand-alone mode
together with the freely available MCR runtime (does not require
MATLAB). Stand-alone mode is supported for 64-bit Windows systems (and
virtual machines VMs) with 4 GB or more RAM memory.

To install, download the installer from link below and run it. The
installer will automatically download and install the correct MCR
runtime engine (>=1 GB file size) and then install the FEATool
executable. Click on the FEATool icon after the installation has been
completed to start the toolbox.<br><br>

[FEATool Multiphysics - Stand-Alone Toolbox Download](https://github.com/precise-simulation/featool-multiphysics/releases/download/1.11.1/FEATool_Multiphysics_webinstaller.exe)

Please note that as the MCR runtime is quite slow to start, the
_FEATool_ GUI may take some minutes to show up after the initial
splash screen has disappeared.


MATLAB Toolbox
--------------

If you have MATLAB, the _FEATool_ toolbox can also be installed directly
from the MATLAB APPS and Add-On Toolbar, the MathWorks File Exchange,
or downloaded directly from the link below.

_FEATool Multiphysics_ has been verified work with 64-bit Windows, Mac
OSX, and Linux operating systems running MATLAB versions 7.9 (R2009b)
and later. Furthermore, a system with 4 GB or more RAM memory is
recommended.

### Installation with MATLAB

- For MATLAB 2012b and later double click on the **[FEATool Multiphysics.mlappinstall](https://github.com/precise-simulation/featool-multiphysics/blob/master/FEATool%20Multiphysics.mlappinstall?raw=true)**
  file, or use the _Get More Apps_ button in the MATLAB _APPS_
  toolbar. Once the app has been installed, a corresponding icon will
  be available in the toolbar to start _FEATool_. (Note that MATLAB
  may not show or give any indication of the app installation progress
  or completion.)

- For MATLAB 2009b-2012a, use the _addpath_ command to add the
  extracted _featool_ program directory to the MATLAB search paths, so
  that the program files can be found by the interpreter (for example
  `addpath C:\featool`). Then simply enter the command `featool` at
  the MATLAB command prompt to start the toolbox GUI and application.

Please note, that using spaces in user and installation directory
paths is not recommended as it can potentially cause issues with
interfaces to external tools, such as geometry engines, grid
generation, and external solvers. Moreover, as all functions are
initially loaded into memory, _FEATool_ may take some time to load and
show the GUI on initial startup.


Tutorials and Examples
----------------------

Pre-defined automated modeling tutorials and examples for various
multi-physics applications can be selected and run from the **File** >
**Model Examples and Tutorials** menu option. Example m-script files
and simulation models are also available in the _examples folder_ of
the _FEATool_ program directory. Moreover, more tutorials and articles
are published on the
[FEATool Technical Articles Blog](https://www.featool.com/articles).


Documentation
-------------

The full
[FEATool Multiphysics Documentation Suite](https://www.featool.com/doc)
is available online and from selecting the corresponding option under
the _Help_ menu in the _FEATool_ GUI.


Basic Use
---------

_FEATool_ and its GUI has been specifically designed to be as easy to
use as possible, and making learning multiphysics simulation by
experimentation possible.

The standard modeling process is divided into six different steps or modes

- **Geometry** - Definition of the geometry to be modeled
- **Grid** - Subdivision of the geometry into smaller cells suitable
  for computation
- **Equation** - Specification of material parameters and coefficients
- **Boundary** - Boundary conditions specify how the model interacts
  with the surrounding environment (outside the geometry)
- **Solve** - Solution and simulation of the defined model problem
- **Post** - Visualization and postprocessing of simulation results

These modes can be accessed by clicking on the corresponding buttons
in left hand side _Mode_ toolbar. The different modes may have
specialized and different _Tools_ available in the corresponding
toolbar. Advanced mode options may also be available in the
corresponding menus.

Basic usage and how to set up and model laminar flow past a cylinder
is explained in the <a href="https://youtu.be/ZnnXl7ryBMI"
target="_blank">linked video tutorial</a> (click on the image to view
the tutorial).

<p align="center"><a href="https://www.youtube.com/watch?v=ZnnXl7ryBMI" target="_blank">
    <img src="https://www.featool.com/images/300-featool-multiphysics-video-tutorial-play.jpg" alt="FEATool Multiphysics Video Tutorial" style="max-width:100%">
</a></p>

The _FEATool_ installation can be tested and validated by running the
simulation model and example test suites by starting _FEATool_ with
either of these commands

    featool testt    % Run tests for GUI tutorials
    featool teste    % Run tests for m-script examples
    featool test     % Run all test suites

Please note that running the test suites may take a significant
amount of time.


External Solvers (_Optional_)
-----------------------------

### OpenFOAM® CFD Solver

The optional OpenFOAM computational fluid dynamics solver integration
makes it easy to perform both laminar and turbulent high performance
CFD simulations directly in MATLAB. OpenFOAM flow simulations often
results in a magnitude or more speedup for instationary simulations
compared to the built-in flow solvers. Additionally, with the
multi-simulation solver integration in _FEATool_ it is possible to
compare and better validate simulation results obtained using both the
built-in, FEniCS, Firedrake, and OpenFOAM solvers.

The OpenFOAM solver binaries are not included with the _FEATool_
distribution and must be installed separately. The FEATool-OpenFOAM
solver integration has been verified with OpenFOAM version 5. For
Windows systems it is recommended to install and use the pre-compiled
blueCFD-core (2017) binaries from
[blueCAPE](http://bluecfd.github.io/Core). For Linux and MacOS systems
the distribution from the
[OpenFOAM Foundation](https://openfoam.org/download) is
recommended. It is necessary that the _simpleFoam_, _pimpleFoam_,
_rhoCentralFoam_, _potentialFoam_, and _collapseEdges_ binaries are
installed and properly set up so they can be called from system script
files (bash scripts on Linux and MacOS and bat/vbs on Windows).


### FEniCS and Firedrake PDE Solvers

[FEniCS](https://fenicsproject.org) and
[Firedrake](https://firedrakeproject.org/) are flexible and
comprehensive finite element analysis (FEA) and partial differential
equation (PDE) modeling and simulation toolkit with Python and C++
interfaces along with many integrated solvers. As both _FEATool_ and
FEniCS discretize equations employing a weak finite element
formulation it is quite straightforward to translate _FEATool_ syntax
and convert it to Python scripts. The FEATool-FEniCS/Firedrake
integration allows for easy conversion, exporting, solving, and
importing _FEATool Multiphysics_ models to FEniCS and Firedrake
directly from the GUI, as well as the MATLAB CLI.

For Ubuntu based Linux and Windows 10 with
[Windows Subsystem for Linux](https://msdn.microsoft.com/en-us/commandline/wsl/about)
systems, FEniCS can be installed by opening a terminal shell and
running the following commands (which automatically also installs
required dependencies such as the Python programming language
interpreter)

    sudo apt-get install software-properties-common
    sudo add-apt-repository ppa:fenics-packages/fenics
    sudo apt-get update
    sudo apt-get install --no-install-recommends fenics

Note that the _FEATool_ distribution does not include a Python
interpreter, FEniCS, or Firedrake itself which also must be installed
separately.  The FEniCS homepage provides instructions
[how to install FEniCS on Linux systems](https://fenicsproject.org/download/)
(note that Docker images are currently not supported by _FEATool_).


External Grid Generators (_Optional_)
-------------------------------------

_FEATool Multiphysics_ comes with built-in support for the external
grid and mesh generators [_Gmsh_](http://gmsh.info/), _Gridgen2D_
(included), and
[_Triangle_](https://www.cs.cmu.edu/~quake/triangle.html), upon
selecting any of these grid generators in the _Grid Settings_ dialog
box, the corresponding binaries will automatically be downloaded and
installed if an internet connection is available. Alternatively, the
mesh generator binaries can downloaded from the
[external grid generators repository](https://github.com/precise-simulation/external-grid-generators)
and/or compiled manually and added to the _FEATool_ installation
directory.

Advantages of using either _Gmsh_, _Gridgen2D_, or _Triangle_ compared
to the built-in grid generation functions is both robustness and mesh
generation speed. Moreover, external grid generators also supports
better and more control with a selection of different mesh generation
algorithms, and specifying the grid size in different geometry
regions, subdomains, as well as on boundaries, allowing for greater
flexibility and better grids tuned for the specific problems and
geometries.

In addition to the external grid generator interfaces, _FEATool_ also
fully supports mesh import and export from the Dolfin/FEniCS (XML),
[GiD](https://www.gidhome.com),
[General Mesh Viewer](http://www.generalmeshviewer.com) (GMV), and
[ParaView](https://www.paraview.org) (VTK/VTP) formats.


License
-------

(C) Copyright 2013-2019 by Precise Simulation Ltd.
All Rights Reserved.

FEATool™ and FEATool Multiphysics™ are trademarks of Precise
Simulation Limited. MATLAB® is a registered trademark of The
MathWorks, Inc.  OPENFOAM® is a registered trade mark of OpenCFD
Limited, producer and distributor of the OpenFOAM® software. All other
trademarks are the property of their respective owners. Precise
Simulation Ltd and its products are not affiliated with, endorsed,
sponsored, or supported by these trademark owners.

The license agreement for using FEATool Multiphysics™ is included with
the distribution and can also be accessed from the _Help_ menu in the
application.

Carefully read the license terms and conditions before installing or
using the programs or documentation. Installing or using the programs
means you have accepted and agree to be bound by the terms and
conditions of this agreement. if you do not accept them, uninstall,
remove and completely delete the programs and documentation.
