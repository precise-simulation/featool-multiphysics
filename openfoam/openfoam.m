%OPENFOAM OpenFOAM CFD solver interface.
%
%   [ U, TLIST, VARS ] = OPENFOAM( FEA, VARARGIN ) Export, solves,
%   and/or imports the solved problem described in the finite
%   element problem struct FEA using the OpenFOAM CFD solver.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}              Description
%       -----------------------------------------------------------------------------------
%       mode        check, export, solve, import Command mode(s) to call (default all)
%       casedir     {tempdir/random}             OpenFOAM case directory
%       control     logical {false}              Show solver control panel.
%       foamdir     default                      OpenFOAM installation directory
%       logfname    default                      OpenFOAM log/output filename
%       fid/logfid  scalar {1/stdout}            Log file/message output file handle
%
%   MODE is a string or cell array of strings selecting action(s) to
%   perform. By default check, export, solve, and import are performed
%   in sequence.
%
%   Returns the solution vector U (n_dof x n_timesteps), corresponding
%   list of time steps TLIST, and additional solution variables in VARS.
%
%   Additional options are passed to the OPENFOAM_DATA,
%   OPENFOAM_EXPORT, OPENFOAM_SOLVE, and OPENFOAM_IMPORT functions.
%
%   Examples:
%
%      1) Laminar steady Hagen-Poiseuille flow in a channel.
%
%      n = 20; rho = 1; miu = 1; uin = 1;
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_rectangle(0, 3, 0, 1)};
%      fea.grid = rectgrid(3*n, 1*n, [0, 3;0, 1]);
%
%      fea = addphys(fea,@navierstokes);
%      fea.phys.ns.eqn.coef{1,end} = {rho};
%      fea.phys.ns.eqn.coef{2,end} = {miu};
%      fea.phys.ns.eqn.coef{5,end} = {uin};
%      fea.phys.ns.bdr.sel(2) = 4;
%      fea.phys.ns.bdr.sel(4) = 2;
%      fea.phys.ns.bdr.coef{2,end}{1,4} = uin;
%
%      fea = parsephys(fea);
%      fea = parseprob(fea);
%
%      fea.sol.u = openfoam(fea);
%
%      subplot(2,1,1)
%      postplot(fea, 'surfexpr', 'p', 'isoexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u', 'v'})
%
%      subplot(2,1,2), hold on, grid on
%      xlabel('Velocity profile at outlet'), ylabel('y')
%      x = 3*ones(1, 100);
%      y = linspace(0, 1, 100);
%      U_ref = 6*uin*(y.*(1-y))./1^2;
%      U = evalexpr('sqrt(u^2+v^2)', [x;y], fea);
%      plot(U_ref, y, 'r--', 'linewidth', 3)
%      plot(U, y, 'b-', 'linewidth', 2.5)
%      legend('Analytic solution', 'Computed solution')
%
%      2) Axisymmetric turbulent flow in a pipe, showing solution convergence curves.
%
%      Re = 1e5; rho = 1; miu = 1/Re; win = 1;
%
%      fea.sdim = {'r', 'z'};
%      fea.geom.objects = {gobj_rectangle(0, .5, 0, 15)};
%      n_lev = 3;
%      nx = 2^(n_lev-1) * 5;
%      ny = 2^(n_lev-1) * 50;
%      px = [.5, .49, .47, .44, .4, .2, 0];
%      px = interp1(linspace(0,0.5,length(px)), px, linspace(0,0.5,nx));
%      fea.grid = rectgrid(px, ny, [0, .5; 0, 15] );
%
%      fea = addphys(fea,{@navierstokes,true});
%      fea.phys.ns.eqn.coef{1,end} = {rho};
%      fea.phys.ns.eqn.coef{2,end} = {miu};
%      fea.phys.ns.eqn.coef{6,end} = {win};
%      fea.phys.ns.bdr.sel(1) = 2;
%      fea.phys.ns.bdr.sel(2) = 1;
%      fea.phys.ns.bdr.sel(3) = 4;
%      fea.phys.ns.bdr.sel(4) = 5;
%      fea.phys.ns.bdr.coef{2,end}{2,1} = win;
%
%      fea = parsephys(fea);
%      fea = parseprob(fea);
%
%      turb.model = 'kEpsilon';
%      turb.inlet = [0.001, 0.00045];
%      turb.wallfcn = 1;
%      fea.sol.u = openfoam(fea, 'turb', turb, 'hax', axes(), 'control', true, 'nproc', 1);
%
%      figure,subplot(1,2,1)
%      postplot(fea, 'surfexpr', 'sqrt(u^2+w^2)', 'isoexpr', 'sqrt(u^2+w^2)', 'arrowexpr', {'u' 'w'})
%      axis([0, .5, 14, 15])
%
%      subplot(1,2,2), hold on, grid on
%      xlabel('Velocity profile at outlet'), ylabel('r')
%      r = linspace(0, 0.5, 100);
%      z = 15*ones(1, 100);
%      U = evalexpr('sqrt(u^2+w^2)', [r;z], fea);
%      plot(U, r, 'b-', 'linewidth', 2.5)
%
%   Further OpenFOAM supported script model examples can be found as
%   EX_NAVIERSTOKES1-4/6-8/10-13/17, EX_COMPRESSIBLEEULER2-6, EX_HEATTRANSFER10.
%
%   See also OPENFOAM_DATA, OPENFOAM_EXPORT, OPENFOAM_SOLVE, OPENFOAM_IMPORT

% Copyright 2013-2024 Precise Simulation, Ltd.
