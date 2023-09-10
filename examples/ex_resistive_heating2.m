function [ fea, out ] = ex_resistive_heating2( varargin )
%EX_RESISTIVE_HEATING2 Resistive heating in a conducting plate.
%
%   [ FEA, OUT ] = EX_RESISTIVE_HEATING2( VARARGIN ) Nonlinear heating
%   in a conductive plate with a time periodic current and temperature
%   dependent resistivity. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       ischeme     scalar 1{2}            Time stepping scheme (1=BE, 2=CN)
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2023 Precise Simulation, Ltd.

cOptDef = { 'ischeme',  2;
            'sfun',     'sflag1';
            'iplot',    1;
            'tol',      [0.03, 0.07, 0.01];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;

% 2D geometry.
p = [-0.05 -0.025; 0 -0.025; 0 -0.005; 0.05 -0.005; 0.05 0.005;
     0 0.005; 0 0.025; -0.05 0.025; -0.05 0.015; -0.01 0.015; -0.01 0.005;
     -0.05 0.005; -0.05 -0.005; -0.01 -0.005; -0.01 -0.015; -0.05 -0.015];
geom.objects{1} = gobj_polygon(p);

% 2D to 3D extrusion.
geom = geom_extrude_face(geom, 'P1', 1, 1e-3, [0,0,1]);

fea.sdim = {'x', 'y', 'z'};
fea.geom.objects = geom.objects(2);

% Mesh generation.
fea.grid = gridgen(fea, 'hmax', 0.0025, 'fid', fid);

% Physics mode for electric potential.
fea = addphys(fea, @conductivemediadc);
fea.phys.dc.sfun = {opt.sfun};
fea.phys.dc.eqn.coef{2,end} = {'s_coef*(1+0.0044*(T-20))'};   % Conductivity.

% Boundary conditions (time dependent input).
fea.phys.dc.bdr.sel(4) = 1;
fea.phys.dc.bdr.coef{2,end}{8} = '100/(0.01*0.001)*sin(2*pi*1*t+(10)/180*pi)';
fea.phys.dc.bdr.coef{2,end}{12} = '100/(0.01*0.001)*sin(2*pi*1*t+(10-120)/180*pi)';
fea.phys.dc.bdr.coef{2,end}{16} = '100/(0.01*0.001)*sin(2*pi*1*t+(10-240)/180*pi)';

% Heat transfer physics mode
fea = addphys(fea, @heattransfer);
fea.phys.ht.sfun = {opt.sfun};
fea.phys.ht.eqn.coef{1,end} = {'rho_coef'};
fea.phys.ht.eqn.coef{2,end} = {'c_coef'};
fea.phys.ht.eqn.coef{3,end} = {'k_coef'};
fea.phys.ht.eqn.coef{7,end} = {'P'};
fea.phys.ht.bdr.sel(:) = 3;
fea.phys.ht.bdr.sel(4) = 1;

% Constants and expressions.
fea.expr = {'s_coef',   '1/1.58e-8';
            'rho_coef', '8900';
            'c_coef',   '385/10';
            'k_coef',   '391';
            'P', '(Vx^2+Vy^2+Vz^2)/(1.58e-8*(1+0.0044*(T-20)))'};

fea = parsephys(fea);
fea = parseprob(fea);

[fea.sol.u,fea.sol.t] = solvetime(fea, 'dt', 0.3, 'tmax', 20, ...
                                  'maxnit', 5, 'ischeme', opt.ischeme, ...
                                  'fid', fid);


% Postprocessing.
if( opt.iplot>0 )
  subplot(2,1,1)
  postplot( fea, 'surfexpr', 'V' )
  title( 'Electric potential, V')

  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'T' )
  title( 'Temperature, T')
end


% Error checking.
[Vmin,Vmax] = minmaxsubd('V', fea);
[Tmin,Tmax] = minmaxsubd('T', fea);

ref = [-6.684e-3, 6.446e-3, 19.233];
out.err = abs(ref - [Vmin,Vmax,Tmax])./abs(ref);
out.pass = all(out.err <= opt.tol) & abs(Tmin)<1e-6;
