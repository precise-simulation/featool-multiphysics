FEATool Multiphysics™ - _Physics Simulation Made Easy_
======================================================

![FEATool Multiphysics Screenshot](https://raw.githubusercontent.com/precise-simulation/featool-multiphysics/master/featool-multiphysics-screenshot.png)

About
-----

[**FEATool Multiphysics**](https://www.featool.com)
(short for Finite Element Analysis Toolbox), is a fully integrated
toolbox for modeling and simulation of coupled physics phenomena,
partial differential equations (PDE), continuum mechanics and
engineering problems.

_FEATool Multiphysics_ aims to provide a truly **easy to use** and
comprehensive **all-in-one** integrated simulation platform for all
kinds of multi-physics analysis. By combining the best of intuitive
**usability**, **extensibility**, and **customization** features,
such as

- Graphical User Interface (GUI)
- Built-in geometry and CAD modeling tools
- Automatic grid and mesh generation
- Postprocessing and visualization
- Fully programmable and scriptable

makes _FEATool Multiphysics_ a suitable simulation and modeling tool
for everyone from students learning mathematical modeling, to
researchers and engineers wishing to explore new ideas in a simple,
fast, and convenient way.


[Features](https://www.featool.com/featool-multiphysics-features)
--------

- Easy to use Graphical User Interface (GUI)
- Built-in and fully integrated geometry and CAD modeling tools, mesh
  generation, multiphysics solvers, post-processing and visualization
- Pre-defined equations and [multi-physics modes](https://www.featool.com/doc/physics.html#phys_modes)
  + [Heat and mass transfer](https://www.featool.com/multiphysics/#heat-and-mass-transfer)
  + [Computational Fluid dynamics (CFD)](https://www.featool.com/cfd-toolbox)
  + [Structural mechanics](https://www.featool.com/multiphysics/#structural-mechanics)
  + [Electromagnetics](https://www.featool.com/multiphysics/#electromagnetics)
  + [Classical PDE](https://www.featool.com/multiphysics/#partial-differential-equations)
- One-click seamless interfaces to external physics solvers
  + [FEniCS (Multiphysics/FEA)](https://www.featool.com/doc/fenics.html)
  + [OpenFOAM® (CFD)](https://www.featool.com/doc/openfoam.html)
  + [SU2 (CFD)](https://www.featool.com/doc/su2.html)
- Full support for [custom and user-defined PDEs](https://www.featool.com/doc/physics.html#phys_ce)
- Equation and expression parser (enter equations and coefficients
  _as writing equations on paper_ without any programming)
- Process, export, and share results and data online with
  [ParaView and Plotly web plots](https://www.featool.com/web-plots)
- Save and export models in
  + Binary file formats
  + MATLAB® m-file script format
  + GUI playback script format
- Fully programmable and scriptable with MATLAB®
  (including support for integration and embedding of
  simulation apps in custom applications and toolboxes)


[System Requirements](https://www.featool.com/doc/quickstart.html#prereq)
-------------------

The _FEATool Multiphysics_ toolbox is a fully integrated simulation
environment for use with MATLAB®. _FEATool_ has been tested and
verified to work with 64-bit Windows, Mac OSX, and Linux operating
systems running MATLAB versions 7.9 (R2009b) and later. Furthermore,
a system with a minimum of 4 GB RAM memory is recommended.


[Installation](https://www.featool.com/doc/quickstart.html#install)
------------

In order to use the _FEATool Multiphysics_ toolbox it must first be
installed on the intended computer system. The toolbox can be
installed from the MATLAB® APPS and Add-On Toolbar, or downloaded
directly from the [Precise Simulation GitHub repository](https://github.com/precise-simulation/featool-multiphysics/releases/latest)
and installed manually.

<p align="center">
  <a href="https://github.com/precise-simulation/featool-multiphysics/raw/master/FEATool%20Multiphysics.mlappinstall" target="_blank"><img src="https://raw.githubusercontent.com/precise-simulation/featool-multiphysics/master/featool-multiphysics-download.png" alt="FEATool Multiphysics Download" style="max-width:50%"></a>
</p>

Please follow the steps below for your matching system to install
_FEATool_ as a MATLAB toolbox. It is recommended to first [uninstall](https://www.featool.com/doc/quickstart.html#uninstall)
previous versions of toolboxes before installing/upgrading to a
newer version. Also note that, as all functionality is loaded into
memory at startup, the toolbox may take some time to fully initialize
and launch the GUI.

### MATLAB 2012b and later

1) First download the latest [FEATool Multiphysics.mlappinstall](https://github.com/precise-simulation/featool-multiphysics/raw/master/FEATool%20Multiphysics.mlappinstall)
   toolbox installation file (if you have downloaded a _zip_ archive
   instead, then extract the _mlappinstall_ file from inside it).

2) Then start MATLAB, press the **APPS** toolbar button,
   and select the **Install App** button.

3) When prompted to choose a toolbox file to install, select the
   **FEATool Multiphysics.mlappinstall** file and press **OK**.

4) Press the **Install** button if prompted to _"Install to My Apps"_.

![FEATool Multiphysics MATLAB Toolbox Installation](https://www.featool.com/doc/featool-multiphysics-toolbox-installation_50.jpg)

Once the toolbox has been installed, an app icon will be available in
the _APPS_ toolbar to start the _FEATool_ GUI. (Note that MATLAB may
not show or give any indication of the toolbox installation progress
or completion.)

### MATLAB 2009b-2012a

1) First download the latest [FEATool Multiphysics zip archive](https://github.com/precise-simulation/featool-multiphysics/archive/master.zip)
   and extract it to a convenient folder on your system.

2) Start MATLAB in this folder, or change to the folder location in
   the MATLAB CLI interpreter with the command
   `cd path_to_featool_program_folder`

3) Run the command `start_featool_gui` from the FEATool program folder
   to start the GUI.

For convenience, one can also use the `addpath path_to_featool_program_folder`
command to permanently add the FEATool program folder to the MATLAB
search paths (one can then start FEATool from any location).


[Tutorials and Examples](https://www.featool.com/doc/quickstart.html#tutorials_and_examples)
----------------------

Pre-defined automated modeling tutorials and examples for various
multi-physics applications can be selected and run from the
          **File** > **Model Examples and Tutorials**
menu option in the GUI.

Example script files and simulation models are also available in the
                       [_examples folder_](https://github.com/precise-simulation/featool-multiphysics/tree/master/examples)
of the _FEATool_ program directory. Moreover, new tutorials and
articles are periodically published on the [FEATool Technical Articles Blog](https://www.featool.com/articles)


[Basic Use](https://www.featool.com/doc/quickstart.html#qs_work)
---------

_FEATool Multiphysics_ and the GUI has been specifically designed to
be as easy to use as possible, and making learning multiphysics
simulation by experimentation possible.

The standard modeling process is divided into six different steps or modes

- **Geometry** - Definition of the geometry to be modeled
- **Grid**     - Subdivision of the geometry into smaller grid cells
                 suitable for computation
- **Equation** - Specification of physics, material parameters, and coefficients
- **Boundary** - Boundary conditions specify how the model interacts
                 with the surrounding environment (outside of the geometry)
- **Solve**    - Solution and simulation of the defined model problem
- **Post**     - Visualization and postprocessing of simulation results

These modes can be accessed by clicking on the corresponding buttons
in left hand side _Mode_ toolbar. Each mode has specialized and
different _Tools_ available in the toolbar that will be activated when
selected. Additional and advanced mode options are also be available
in the mode menus.

Basic usage and how to set up and model coupled fluid flow and
temperature in a heat exchanger is explained in the
[linked video tutorial](https://youtu.be/TBfVWgYbGTw)
(click on the image below to start the tutorial).

<p align="center">
  <a href="https://www.youtube.com/watch?v=TBfVWgYbGTw" target="_blank">
    <img src="https://img.youtube.com/vi/TBfVWgYbGTw/0.jpg"
         alt="FEATool Multiphysics Video Tutorial" style="max-width:100%">
  </a>
</p>


Documentation
-------------

The full [FEATool Multiphysics Documentation Suite](https://www.featool.com/doc)
is available online, and by selecting the corresponding
option in the _Help_ menu of the _FEATool_ GUI.


License
-------

(C) Copyright 2013-2021 by Precise Simulation Ltd.
All Rights Reserved.

FEATool™ and FEATool Multiphysics™ are trademarks of Precise
Simulation Limited. MATLAB® is a registered trademark of The
MathWorks, Inc. OPENFOAM® is a registered trade mark of OpenCFD
Limited, producer and distributor of the OpenFOAM® software.
All other trademarks are the property of their respective
owners. Precise Simulation Limited and its products are not affiliated
with, endorsed, sponsored, or supported by these trademark owners.

The license agreement for using FEATool Multiphysics™ is included with
the distribution and can also be viewed by selecting
_About FEATool..._ > _License Agreement_ from the _Help_ menu
in the application.

Carefully read the license terms and conditions before installing or
using the programs or documentation. Installing or using the programs
means you have accepted and agree to be bound by the terms and
conditions of this agreement. if you do not accept them, uninstall,
remove and completely delete the programs and documentation.
