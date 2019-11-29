function featool()
% FEATool Multiphysics™ - _Physics Simulation Made Easy_
% ======================================================
% 
% ![FEATool Multiphysics Screenshot](featool-multiphysics-screenshot.png)
% 
% 
% About
% -----
% 
% **FEATool Multiphysics™** (short for Finite Element Analysis Toolbox),
% is a toolbox for modeling and simulation of coupled physical
% phenomena, [partial differential equations](https://en.wikipedia.org/wiki/Partial_differential_equation)
% (PDE), continuum mechanics and engineering problems with the
% [finite element method](https://en.wikipedia.org/wiki/Finite_element_method)
% (FEM).
% 
% _FEATool_ aims to provide an easy to use and comprehensive all-in-one
% integrated simulation package for all types of finite element and
% multiphysics analyses, and combine the best of ease of use, powerful
% functionality, and extensibility for both beginners and advanced
% users. Features such as an intuitive and easy to learn graphical user
% interface (GUI), a complete library of grid generation, and
% postprocessing functions, as well as command line interface (CLI)
% programming, and interactive and interpreted programming and scripting
% capabilities, makes _FEATool_ suitable for everyone from students
% learning mathematical modeling, to professionals and engineers wishing
% to explore new ideas in a simple, fast, and convenient way.
% 
% 
% Key Features
% ------------
% 
% - Easy to use Graphical User Interface (GUI)
% - Fully integrated and built-in pre-processing, CAD tools, automatic
%   mesh generation, solvers, post-processing, and visualization
%   functionality
% - Pre-defined physics modes and equations for heat and mass transfer,
%   fluid dynamics, structural mechanics, electromagnetics, and
%   classical PDE
% - Support for custom user-defined equations and PDEs
% - Built-in expression parser (enter equations and coefficients _just
%   as writing on paper_ without any programming)
% - Easy to define fully coupled multiphysics problems
% - One-click interfaces for the OpenFOAM®, FEniCS, and Firedrake
%   external solvers
% - Save and export models in binary file format or editable script
%   files (every GUI operation is recorded to a corresponding MATLAB
%   function call)
% - Fully scriptable and programmable with the MATLAB scripting language
%   (including support for integration and embedding in custom
%   applications and toolboxes)
% 
% 
% System Requirements
% -------------------
% 
% The _FEATool Multiphysics_ Toolbox can either be used stand-alone as a
% fully integrated simulation GUI, or together with MATLAB enabling
% command line (CLI) use and scripting. In order to use the _FEATool
% Multiphysics_ toolbox, the program files must first be installed on the
% intended computer system. Please follow the procedure below for the
% corresponding stand-alone or MATLAB toolbox installation modes. Note
% that it is recommended to first uninstall older versions of _FEATool_
% before updating to a newer version.
% 
% 
% Stand-alone Toolbox
% -------------------
% 
% FEATool Multiphysics has been compiled to run in stand-alone mode
% together with the freely available MCR runtime (does not require
% MATLAB). Stand-alone mode is supported for 64-bit Windows systems (and
% virtual machines VMs) with 4 GB or more RAM memory.
% 
% To install, download the installer from link below and run it. The
% installer will automatically download and install the correct MCR
% runtime engine (>=1 GB file size) and then install the FEATool
% executable. Click on the FEATool icon after the installation has been
% completed to start the toolbox.<br><br>
% 
% [![FEATool Multiphysics - Stand-Alone Toolbox Download](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/multiphysics/files/1.11.1/FEATool_Multiphysics_webinstaller.exe/download)<br>
% 
% Please note that as the MCR runtime is quite slow to start, the
% _FEATool_ GUI may take some minutes to show up after the initial
% splash screen has disappeared.
% 
% 
% MATLAB Toolbox
% --------------
% 
% If you have MATLAB, the _FEATool_ toolbox can also be installed directly
% from the MATLAB APPS and Add-On Toolbar, the MathWorks File Exchange,
% or downloaded directly from the link below.
% 
% _FEATool Multiphysics_ has been verified work with 64-bit Windows, Mac
% OSX, and Linux operating systems running MATLAB versions 7.9 (R2009b)
% and later. Furthermore, a system with 4 GB or more RAM memory is
% recommended.
% 
% ### Installation with MATLAB
% 
% - For MATLAB 2012b and later double click on the **FEATool Multiphysics.mlappinstall**
%   file, or use the _Get More Apps_ button in the MATLAB _APPS_
%   toolbar. Once the app has been installed, a corresponding icon will
%   be available in the toolbar to start _FEATool_. (Note that MATLAB
%   may not show or give any indication of the app installation progress
%   or completion.)
% 
% - For MATLAB 2009b-2012a, use the _addpath_ command to add the
%   extracted _featool_ program directory to the MATLAB search paths, so
%   that the program files can be found by the interpreter (for example
%   `addpath C:\featool`). Then simply enter the command `featool` at
%   the MATLAB command prompt to start the toolbox GUI and application.
% 
% Please note, that using spaces in user and installation directory
% paths is not recommended as it can potentially cause issues with
% interfaces to external tools, such as geometry engines, grid
% generation, and external solvers. Moreover, as all functions are
% initially loaded into memory, _FEATool_ may take some time to load and
% show the GUI on initial startup.
% 
% 
% Tutorials and Examples
% ----------------------
% 
% Pre-defined automated modeling tutorials and examples for various
% multi-physics applications can be selected and run from the **File** >
% **Model Examples and Tutorials** menu option. Example m-script files
% and simulation models are also available in the _examples folder_ of
% the _FEATool_ program directory. Moreover, more tutorials and articles
% are published on the
% [FEATool Technical Articles Blog](https://www.featool.com/articles).
% 
% 
% Documentation
% -------------
% 
% The full
% [FEATool Multiphysics Documentation Suite](https://www.featool.com/doc)
% is available online and from selecting the corresponding option under
% the _Help_ menu in the _FEATool_ GUI.
% 
% 
% Basic Use
% ---------
% 
% _FEATool_ and its GUI has been specifically designed to be as easy to
% use as possible, and making learning multiphysics simulation by
% experimentation possible.
% 
% The standard modeling process is divided into six different steps or modes
% 
% - **Geometry** - Definition of the geometry to be modeled
% - **Grid** - Subdivision of the geometry into smaller cells suitable
%   for computation
% - **Equation** - Specification of material parameters and coefficients
% - **Boundary** - Boundary conditions specify how the model interacts
%   with the surrounding environment (outside the geometry)
% - **Solve** - Solution and simulation of the defined model problem
% - **Post** - Visualization and postprocessing of simulation results
% 
% These modes can be accessed by clicking on the corresponding buttons
% in left hand side _Mode_ toolbar. The different modes may have
% specialized and different _Tools_ available in the corresponding
% toolbar. Advanced mode options may also be available in the
% corresponding menus.
% 
% Basic usage and how to set up and model laminar flow past a cylinder
% is explained in the <a href="https://youtu.be/ZnnXl7ryBMI"
% target="_blank">linked video tutorial</a> (click on the image to view
% the tutorial).
% 
% <p align="center"><a href="https://www.youtube.com/watch?v=ZnnXl7ryBMI" target="_blank">
%     <img src="https://www.featool.com/images/300-featool-multiphysics-video-tutorial-play.jpg" alt="FEATool Multiphysics Video Tutorial" style="max-width:100%">
% </a></p>
% 
% The _FEATool_ installation can be tested and validated by running the
% simulation model and example test suites by starting _FEATool_ with
% either of these commands
% 
%     featool testt    % Run tests for GUI tutorials
%     featool teste    % Run tests for m-script examples
%     featool test     % Run all test suites
% 
% Please note that running the test suites may take a significant
% amount of time.
% 
% 
% External Solvers (_Optional_)
% -----------------------------
% 
% ### OpenFOAM® CFD Solver
% 
% The optional OpenFOAM computational fluid dynamics solver integration
% makes it easy to perform both laminar and turbulent high performance
% CFD simulations directly in MATLAB. OpenFOAM flow simulations often
% results in a magnitude or more speedup for instationary simulations
% compared to the built-in flow solvers. Additionally, with the
% multi-simulation solver integration in _FEATool_ it is possible to
% compare and better validate simulation results obtained using both the
% built-in, FEniCS, Firedrake, and OpenFOAM solvers.
% 
% The OpenFOAM solver binaries are not included with the _FEATool_
% distribution and must be installed separately. The FEATool-OpenFOAM
% solver integration has been verified with OpenFOAM version 5. For
% Windows systems it is recommended to install and use the pre-compiled
% blueCFD-core (2017) binaries from
% [blueCAPE](http://bluecfd.github.io/Core). For Linux and MacOS systems
% the distribution from the
% [OpenFOAM Foundation](https://openfoam.org/download) is
% recommended. It is necessary that the _simpleFoam_, _pimpleFoam_,
% _rhoCentralFoam_, _potentialFoam_, and _collapseEdges_ binaries are
% installed and properly set up so they can be called from system script
% files (bash scripts on Linux and MacOS and bat/vbs on Windows).
% 
% 
% ### FEniCS and Firedrake PDE Solvers
% 
% [FEniCS](https://fenicsproject.org) and
% [Firedrake](https://firedrakeproject.org/) are flexible and
% comprehensive finite element analysis (FEA) and partial differential
% equation (PDE) modeling and simulation toolkit with Python and C++
% interfaces along with many integrated solvers. As both _FEATool_ and
% FEniCS discretize equations employing a weak finite element
% formulation it is quite straightforward to translate _FEATool_ syntax
% and convert it to Python scripts. The FEATool-FEniCS/Firedrake
% integration allows for easy conversion, exporting, solving, and
% importing _FEATool Multiphysics_ models to FEniCS and Firedrake
% directly from the GUI, as well as the MATLAB CLI.
% 
% For Ubuntu based Linux and Windows 10 with
% [Windows Subsystem for Linux](https://msdn.microsoft.com/en-us/commandline/wsl/about)
% systems, FEniCS can be installed by opening a terminal shell and
% running the following commands (which automatically also installs
% required dependencies such as the Python programming language
% interpreter)
% 
%     sudo apt-get install software-properties-common
%     sudo add-apt-repository ppa:fenics-packages/fenics
%     sudo apt-get update
%     sudo apt-get install --no-install-recommends fenics
% 
% Note that the _FEATool_ distribution does not include a Python
% interpreter, FEniCS, or Firedrake itself which also must be installed
% separately.  The FEniCS homepage provides instructions
% [how to install FEniCS on Linux systems](https://fenicsproject.org/download/)
% (note that Docker images are currently not supported by _FEATool_).
% 
% 
% External Grid Generators (_Optional_)
% -------------------------------------
% 
% _FEATool Multiphysics_ comes with built-in support for the external
% grid and mesh generators [_Gmsh_](http://gmsh.info/), _Gridgen2D_
% (included), and
% [_Triangle_](https://www.cs.cmu.edu/~quake/triangle.html), upon
% selecting any of these grid generators in the _Grid Settings_ dialog
% box, the corresponding binaries will automatically be downloaded and
% installed if an internet connection is available. Alternatively, the
% mesh generator binaries can downloaded from the
% [external grid generators repository](https://github.com/precise-simulation/external-grid-generators)
% and/or compiled manually and added to the _FEATool_ installation
% directory.
% 
% Advantages of using either _Gmsh_, _Gridgen2D_, or _Triangle_ compared
% to the built-in grid generation functions is both robustness and mesh
% generation speed. Moreover, external grid generators also supports
% better and more control with a selection of different mesh generation
% algorithms, and specifying the grid size in different geometry
% regions, subdomains, as well as on boundaries, allowing for greater
% flexibility and better grids tuned for the specific problems and
% geometries.
% 
% In addition to the external grid generator interfaces, _FEATool_ also
% fully supports mesh import and export from the Dolfin/FEniCS (XML),
% [GiD](https://www.gidhome.com),
% [General Mesh Viewer](http://www.generalmeshviewer.com) (GMV), and
% [ParaView](https://www.paraview.org) (VTK/VTP) formats.
% 
% 

% License
% -------
% 
% (C) Copyright 2013-2019 by Precise Simulation Ltd.
% All Rights Reserved.
% 
% FEATool™ and FEATool Multiphysics™ are trademarks of Precise
% Simulation Limited. MATLAB® is a registered trademark of The
% MathWorks, Inc.  OPENFOAM® is a registered trade mark of OpenCFD
% Limited, producer and distributor of the OpenFOAM® software. All other
% trademarks are the property of their respective owners. Precise
% Simulation Ltd and its products are not affiliated with, endorsed,
% sponsored, or supported by these trademark owners.
% 
% The license agreement for using FEATool Multiphysics™ is included with
% the distribution and can also be accessed from the _Help_ menu in the
% application.
% 
% Carefully read the license terms and conditions before installing or
% using the programs or documentation. Installing or using the programs
% means you have accepted and agree to be bound by the terms and
% conditions of this agreement. if you do not accept them, uninstall,
% remove and completely delete the programs and documentation.

% Precise Simulation Limited Software License Agreement
% 
% CAREFULLY READ THE FOLLOWING TERMS AND CONDITIONS ("TERMS") BEFORE
% INSTALLING OR USING THE PROGRAMS OR DOCUMENTATION. INSTALLING OR USING
% THE PROGRAMS MEANS YOU HAVE ACCEPTED AND AGREE TO BE BOUND BY THE
% TERMS AND CONDITIONS OF THIS AGREEMENT. IF YOU DO NOT ACCEPT THEM,
% UNINSTALL, REMOVE AND COMPLETELY DELETE THE PROGRAMS AND
% DOCUMENTATION.
% 
% 1. Preamble: This Agreement governs the relationship between the
% Licensee ("you", "your") and Licensor Precise Simulation Limited
% ("we", "us", "ours"), a duly registered company whose registered place
% of business is Suite 601, 6/F, Tai Tung Building, 8 Fleming Road, Wan
% Chai, Hong Kong. This Agreement sets the terms, rights, restrictions
% and obligations on using FEATool ("Software", "Program(s)") and
% documentation ("Documentation") created and owned by Licensor, as
% detailed herein.
% 
% 2. License Grant: Licensor hereby grants Licensee a Non-assignable &
% Non-transferable, Non-exclusive license to run and use the Program,
% without the rights to create derivative works, all with accordance
% with the terms set forth and other legal restrictions set forth in 3rd
% party software used while running Software.
% 
% 2.1 Programs: You may license a specified single installation license
% ("SUL"), multi-user/floating network license ("MUL"), or ("CKL") class
% kit license under this Agreement, and your license rights are for the
% number of installations and users set forth on the purchase order,
% agreement, or issued invoice. A free limited and restricted license
% ("FREE/TRIAL") is granted for personal, non-commercial use for
% evaluation purposes.
% 
% a. the FREE/TRIAL license option is restricted to personal, trial, and
% non-commercial use allowing for a single installation and concurrent
% use of the Program. You may NOT use the Program with a FREE/TRIAL
% license for any commercial, or production use, i.e., you may only use
% the Program for experimental, personal, and trial use (to test the
% Program). Specifically, the restrictions of the FREE/TRIAL license
% Program and Software may not be circumvented in any way without
% Payment for an upgraded license.
% 
% b. the specified single installation license SUL must be installed on
% a specified computer system and its use is limited to a single
% concurrent instance. To change system a system transfer fee may be
% required.
% 
% c. the multi-use license option MUL may be installed on a single
% networked system or server, or several systems and run concurrently
% the number of instances specified in the purchase order, agreement, or
% issued invoice.
% 
% d. academic granting institutions with the class kit license CKL
% option may install and use the Software in a computer lab/systems
% belonging to the institute/institution and run concurrently the number
% of instances specified in the purchase order, agreement, or issued
% invoice.
% 
% e. regardless of which license you have, you shall use the Programs
% only for your internal operations. For the purposes of this Agreement,
% "internal operations" means use of the Programs by your employees or
% those of your subsidiaries or parent company and for the performance
% of consulting or research for third parties who engage you as an
% employee or independent contractor. You also shall not disclose any
% characteristics or technical capabilities of the Programs to any third
% party without our prior written authorization.
% 
% 2.2 Delivery: We may deliver the Programs and Documentation to you in
% archival form over the Internet with a passcode or license key which
% specifies the licensed Programs. You shall be responsible for all use
% of your passcode, authorized or not, and you shall not disclose the
% archive passcode or allow it to be used except for installation of the
% Programs.
% 
% 2.3 Ownership: All right, title and interest in and to the licensed
% Program(s), including without limitation, trade secrets and
% copyrights, are, and shall at all times remain, the exclusive property
% of us and you shall have no right, therein, except the expressly
% limited license rights granted herein.
% 
% 2.4. Non Assignable & Non-Transferable: Licensee may not assign or
% transfer his rights and duties under this license.
% 
% 2.5. The Software and Documentation are for your personal use and/or
% internal business operations and are not for resale or other transfer
% or disposition to any other person or entity. In addition, you
% specifically agree not to:
% 
% a. reverse engineer, decompile, disassemble, translate, modify, alter
% or otherwise change the Licensor's Software or any part thereof;
% 
% b. attempt to derive the source code, design or structure of the
% Licensor's Software;
% 
% c. sell, rent, lease, distribute, assign, sub-license, convey,
% transfer, pledge as security or otherwise encumber or transfer
% (including by loan or gift) the rights and licenses granted hereunder;
% 
% d. copy or reproduce any part of the Software or Documentation other
% than as allowed under this Agreement;
% 
% e. use the Software or Documentation in any manner that violates any
% statute, law, rule, regulation, directive, guideline, bylaw whether
% presently in force or may be implemented by state or local
% authorities.
% 
% 3. Term & Termination: The Term of this license shall be until
% terminated, or until specified by issued purchase order, agreement, or
% issued invoice. Licensor may terminate this Agreement, including
% Licensee's license in the case where Licensee:
% 
% a. became insolvent or otherwise entered into any liquidation process; or
% 
% b. Licensee was in breach of any of this license's terms and
% conditions and such breach was not cured, immediately upon
% notification; or
% 
% c. Licensee otherwise entered into any arrangement which caused
% Licensor to be unable to enforce his rights under this License.
% 
% 4. Payment: In consideration of the License granted under clause 2,
% Licensee shall pay Licensor a fee which Licensor may deem
% adequate. Failure to perform payment shall construe as material breach
% of this Agreement. You shall be liable for any taxes (except those on
% our net income) due in connection with this Agreement.
% 
% 4.1 No purchase order or any other standardized business form issued
% by you, and even if such purchase order or other standardized business
% form provides that it takes precedence over any other agreement
% between the parties, shall be effective to contradict, modify, add to
% or delete from the terms of this Agreement in any manner
% whatsoever. Any acknowledgment, in any form, of any such purchase
% order or standardized business form is not recognized as a subsequent
% writing and will not act as acceptance of such terms.
% 
% 5. Upgrades, Updates and Fixes: Licensor may provide Licensee, from
% time to time, with Upgrades, Updates or Fixes, as detailed herein and
% according to his sole discretion. Licensee hereby warrants to keep The
% Software up-to-date and install all relevant updates and fixes, and
% may, at his sole discretion, purchase upgrades, according to the rates
% set by Licensor. Licensor shall provide any update or Fix free of
% charge; however, nothing in this Agreement shall require Licensor to
% provide Updates or Fixes.
% 
% 6. Support: The Software is provided under an AS-IS basis and without
% any support, updates or maintenance. Nothing in this Agreement shall
% require Licensor to provide Licensee with support or fixes to any bug,
% failure, mis-performance or other defect in The Software.
% 
% 7. Liability: To the extent permitted under Law, The Software is
% provided under an AS-IS basis. Licensor shall never, and without any
% limit, be liable for any damage, cost, expense or any other payment
% incurred by Licensee as a result of Software's actions, failure, bugs
% and/or any other interaction between The Software and Licensee's
% end-equipment, computers, other software or any 3rd party,
% end-equipment, computer or services.  Moreover, Licensor shall never
% be liable for any defect in source code written by Licensee when
% relying on The Software or using The Software's source code.
% 
% 8. Warranty: The Software is provided without any warranty; Licensor
% hereby disclaims any warranty that The Software shall be error free,
% without defects or code which may cause damage to Licensee's computers
% or to Licensee, and that Software shall be functional. Licensee shall
% be solely liable to any damage, defect or loss incurred as a result of
% operating software and undertake the risks contained in running The
% Software on License's Computer System(s) and Server(s).
% 
% 8.1 Prior Inspection: Licensee hereby states that he inspected The
% Software thoroughly and found it satisfactory and adequate to his
% needs, that it does not interfere with his regular operation and that
% it does meet the standards and scope of his computer systems and
% architecture. Licensee found that The Software interacts with his
% development, website and server environment and that it does not
% infringe any of End User License Agreement of any software Licensee
% may use in performing his services. Licensee hereby waives any claims
% regarding The Software's incompatibility, performance, results and
% features, and warrants that he inspected the The Software.
% 
% 9. No Refunds: Licensee warrants that he inspected The Software
% according to clause 8.1 and that it is adequate to his
% needs. Accordingly in the case of NON-FREE licenses, as The Software
% is intangible goods, Licensee shall not be, ever, entitled to any
% refund, rebate, compensation or restitution for any reason whatsoever,
% even if The Software contains material flaws.
% 
% 10. Technical Information. You agree that We may collect or process
% technical and related information arising from Your use of the
% Software which may include but may not be limited to internet protocol
% address, hardware identification, operating system, application
% software, peripheral hardware, and non-personally identifiable
% Software usage statistics to facilitate the provisioning of Updates,
% Support, invoicing or online services, identify trends and bugs,
% collect activation information, usage statistics and track other data
% related to Your use of the Software.
% 
% 11. Indemnification: Licensee hereby warrants to hold Licensor
% harmless and indemnify Licensor for any lawsuit brought against it in
% regards to Licensee's use of The Software in means that violate,
% breach or otherwise circumvent this license, Licensor's intellectual
% property rights or Licensor's title in The Software. Licensor shall
% promptly notify Licensee in case of such legal action and request
% Licensee's consent prior to any settlement in relation to such lawsuit
% or claim.
% 
% 12. Governing Law, Jurisdiction: Licensee hereby agrees not to
% initiate class-action lawsuits against Licensor in relation to this
% license and to compensate Licensor for any legal fees, cost or
% attorney fees should any claim brought by Licensee against Licensor be
% denied, in part or in full.
% 
% 13. Revised Terms of Use: We may revise the terms of use of the
% Programs from time to time. Revisions are effective upon receipt of
% notice from us.
