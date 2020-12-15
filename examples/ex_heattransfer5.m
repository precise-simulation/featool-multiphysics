function [ fea, out ] = ex_heattransfer5( varargin )
%EX_HEATTRANSFER5 2D Transient cooling and shrink fitting example.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER5( VARARGIN ) This example models transient cooling
%   for shrink fitting a two part assembly. A tungsten rod heated to 84 C is inserted
%   into a -10 C chilled steel frame part. The time when the maximum temperature has
%   cooled to 70 C should be determined. The assembly is cooled due to convection
%   through a surrounding medium kept at 17 C and a constant heat transfer coefficient
%   of 50 W/(m^2 K) [1].
%
%   Reference:
%
%      [1] Krysl P. A Pragmatic Introduction to the Finite Element Method
%          for Thermal and Stress Analysis. Pressure Cooker Press, USA, 2005.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.0025}        Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2020 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.005;
            'sfun',     'sflag1';
            'solver',   '';
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry definition.
r1 = gobj_rectangle( 0, 0.11, 0, 0.12,  'R1' );
c1 = gobj_circle( [ 0.065 0 ],   0.015, 'C1' );
c2 = gobj_circle( [ 0.11 0.12 ], 0.035, 'C2' );
c3 = gobj_circle( [ 0 0.06 ],    0.025, 'C3' );

r2 = gobj_rectangle( 0.065, 0.16, 0.05, 0.07, 'R2' );
c4 = gobj_circle( [ 0.065 0.06 ], 0.01, 'C4' );

fea.geom.objects = { r1 c1 c2 c3 r2 c4 };
fea = geom_apply_formula( fea, 'R1-C1-C2-C3' );
fea = geom_apply_formula( fea, 'R2+C4' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );


% Problem definition.
fea.sdim  = { 'x', 'y' };             % Space coordinate names.
fea = addphys( fea, @heattransfer );  % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };      % Set shape function.

% Equation coefficients.
rho_steel    =  7500;
rho_tungsten = 19000;
cp_steel     =   470;
cp_tungsten  =   134;
k_steel      =    44;
k_tungsten   =   163;
fea.phys.ht.eqn.coef{1,end} = { rho_steel rho_tungsten rho_tungsten };   % Density
fea.phys.ht.eqn.coef{2,end} = { cp_steel  cp_tungsten  cp_tungsten  };   % Heat capcity.
fea.phys.ht.eqn.coef{3,end} = { k_steel   k_tungsten   k_tungsten   };   % Thermal conductivity.
fea.phys.ht.eqn.coef{7,end} = { -10 84 84 };                             % Initial temperature.

% Boundary conditions.
n_bdr = max(fea.grid.b(3,:));
fea.phys.ht.bdr.sel(fea.phys.ht.bdr.sel>0) = 4;
h_coef = 50;
for i_bdr=1:n_bdr
  fea.phys.ht.bdr.coef{4,end}{i_bdr}{2} = h_coef;
  fea.phys.ht.bdr.coef{4,end}{i_bdr}{3} = 17;
end


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', 0.25, 'tmax', 16, 'ischeme', 2 );
  tlist = fea.sol.t;
else
  [fea.sol.u, tlist] = solvetime( fea, 'fid', opt.fid, 'tmax', 16, 'tstep', 0.25, 'tstop', 0, ...
                                  'init', {'T0_ht'}, 'maxnit', 2 );
end

% Check when max temperature < 70.
for i=1:numel(tlist)
  T_min(i) = min(fea.sol.u(:,i));
  T_max(i) = max(fea.sol.u(:,i));
end
ix = find(T_max<70);
i1 = ix(1);
i2 = i1 - 1;
s  = ( T_max(i2) - 70 )/( T_max(i2) - T_max(i1) );
t_70 = tlist(i2) + s*( tlist(i1) - tlist(i2) );
u_70 = fea.sol.u(:,i2) + s*( fea.sol.u(:,i1) - fea.sol.u(:,i2) );


% Postprocessing.
if( opt.iplot>0 )
  figure
  fea_plot = fea;
  fea_plot.sol.u = u_70;
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'T', 'isoexpr', 'T', 'isolev', 20 )
  title(['Temperature distribution at t=',num2str(t_70)])

  subplot(1,2,2)
  plot( tlist, T_min, 'b-' )
  hold on
  plot( tlist, T_max, 'r-' )
  grid on
  title('Maximum and minimum temperatures')
  ylabel('Temperature [C]')
  xlabel('Time [s]')
end


% Error checking.
out.T_min = T_min;
out.T_max = T_max;
out.T_70  = t_70;
out.pass = abs(t_70-13)<3.5;


if( nargout==0 )
  clear fea out
end
