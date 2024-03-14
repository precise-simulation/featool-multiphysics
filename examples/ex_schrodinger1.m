function [ fea, out ] = ex_schrodinger1( varargin )
%EX_SCHRODINGER1 Schrodinger equation for the Hydrogen atom
%
%   [ FEA, OUT ] = EX_SCHRODINGER1( VARARGIN ) Computation of energy
%   levels and electron orbits for a 1-particle Hydrogen atom system
%   using the Schrodinger equation in cylindrical coordiantes. Accepts
%   the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.1}           Grid cell size (in meter scale)
%       sfun        string {sflag2}        Finite element shape function
%       iphys       scalar 0/{1}           Use physics mode to define problem    (=1)
%                                          or directly define fea.eqn/bdr fields (=0)
%                                          or use core assembly functions        (<0)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2024 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.05;
            'sfun',     'sflag2';
            'iplot',    true;
            'tol',      0.02;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition (in meters).
fea.sdim = {'r', 'z'};
fea.geom.objects = { gobj_circle([0,0], 3, 'C1'), ...
                     gobj_circle([0,0], 0.5, 'C2'), ...
                     gobj_rectangle(-3, 0, -3, 3, 'R1'), ...
                     gobj_rectangle(-3, 0, -3, 3, 'R2')};
fea = geom_apply_formula( fea, 'C1-R1' );
fea = geom_apply_formula( fea, 'C2-R2' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );

% Scale geometry and grid to nanometer scale.
s = 1e-9;
fea = geom_apply_transformation( fea, fea.geom.tags, 0, s );
fea.grid.p = fea.grid.p * s;


% Constants and expressions.
me = 9.10939e-31;     % Electron mass.
ec = 1.6021773e-19;   % Electron charge.
hp = 6.626076e-34;    % Planck's constant.
e0 = 8.8541878e-12;   % Electric capacitivity in vacuum.
m  = 0;               % Ground state energy.
fea.expr = { 'me', me;
             'e',  ec;
             'hp', hp;
             'e0', e0;
             'C',  '8*me*pi^2/hp^2';
             'm',  m};


% Problem definition.
fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
fea.phys.poi.sfun = { opt.sfun };     % Set FEM shape function.

fea.phys.poi.eqn.seqn = ...
'C*r*u'' - r*(ur_r+uz_z) = 1+(C*r*e^2/(4*pi*e0*sqrt(r^2+z^2))-m^2/r)*u_t';
% fea.phys.poi.eqn.coef{1,end} = { 'C*r' };   % Set time coefficient.
% fea.phys.poi.eqn.coef{2,end} = { 'r' };     % Set diffusion coefficient.
% fea.phys.poi.eqn.coef{3,end} = { '1+(C*r*e^2/(4*pi*e0*sqrt(r^2+z^2))-m^2/r)*u' };       % Define source term coefficient.


% Define boundary conditions (Homogenous Neumann on symmetry axis, and Dirichlet on sphere boundary.)
fea.phys.poi.bdr.sel = [1 1 2 2 2];


% Parse and solve problem.
fea = parsephys(fea);
fea = parseprob(fea);
[fea.sol.u, fea.sol.l] = solveeig( fea, 'sigma', -2.18e-18, ...
                                   'neigs', 10, 'fid', fid );


% Error checking.
lam = [-2.18e-18; -5.45e-19; -5.45e-19; -2.422e-19; -2.422e-19];
out.err = abs((lam - fea.sol.l(1:5))./lam);
out.pass = all(out.err < opt.tol);


% Postprocessing.
if( opt.iplot>0 )
  isol = 5;
  postplot( fea, 'surfexpr', 'u', 'solnum', isol )
  title(sprintf('u(lambda(%d)) = %g', isol, fea.sol.l(isol)))
end

if( nargout==0 )
  clear fea out
end
