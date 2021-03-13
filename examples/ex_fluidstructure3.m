function [ fea, out ] = ex_fluidstructure3( varargin )
%EX_FLUIDSTRUCTURE3 Fluid-structure interaction - elastic beam.
%
%   [ FEA, OUT ] = EX_FLUIDSTRUCTURE3( VARARGIN ) Example for
%   fluid-structure interaction flow around an elastic beam.
%
%   Reference:
%
%   [1] Wall W., Ramm E. Fluid-Interaktion mit stabilisierten Finiten
%   Elementen. Phd Thesis, Institut fur  Baustatik, Universitat
%   Stuttgart, 1999.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry.
fea.sdim = { 'x', 'y' };
gobj1 = gobj_rectangle( 0, 19.5, 0, 12, 'R1' );
gobj2 = gobj_rectangle( 3.5, 13.5, 3.5, 8.5, 'R2' );
gobj3 = gobj_rectangle( 4.5, 5.5, 5.5, 6.5, 'R3' );
fea.geom.objects = { gobj1, gobj2, gobj3 };
fea.geom = copy_geometry_object( 'R3', fea.geom );
gobj = gobj_rectangle( 5.5, 9, 6-0.03, 6+0.03, 'R5' );
fea  = geom_add_gobj( fea, gobj );
fea.geom = geom_apply_formula( fea.geom, 'R1-R3' );
fea.geom = geom_apply_formula( fea.geom, 'R2-R4' );


% Grid generation.
hmaxb = 0.75*ones(1,17);
hmaxb([1:6,11:13]) = 0.015;
hmaxb([14:17]) = 0.2;
fea.grid = gridgen( fea, 'hmaxb', hmaxb, 'gridgen', 'gmsh', 'fid', opt.fid );


% Equation settings.
fea = addphys( fea, @fluidstructure );
fea.phys.fsi.eqn.coef{1,end} = { '0.1', '1.18*1e-3', '1.18*1e-3' };   % Density.
fea.phys.fsi.eqn.coef{2,end} = { '1', '1.82*1e-4', '1.82*1e-4' };     % Viscosity.
fea.phys.fsi.eqn.coef{3,end} = { '0.25', '0.3', '0.3' };              % Poisson's ratio.
fea.phys.fsi.eqn.coef{4,end} = { '2.5e6', '1', '1' };                 % Modulus of elasticity.
fea.phys.fsi.prop.active = [ 0, 1, 1; 0, 1, 1; 0, 1, 1; 1, 0, 0; 1, 0, 0 ];


% Boundary settings.
fea.phys.fsi.bdr.sel = [ 6, 1, 1, 1, 1, 1, 5, 3, 5, 2, -2, -2, -2, -1, -1, -1, -1 ];
fea.phys.fsi.bdr.coef{2,end}{1,10} = 51.3;


% Solver.
fea = parsephys(fea);
fea = parseprob(fea);

[fea.sol.u,fea.sol.t,fea.sol.grid.p] = fsisolve( fea, 'tstep', 0.01, 'tmax', 1, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot(fea,'surfexpr','p')
end


% Error checking.
v_min = inf; v_max = -inf;
for i=80:size(fea.sol.u,2)
  [vm,vx] = minmaxsubd( 'dy', fea, [], [], [], i );
  if( vm<v_min )
    v_min = vm;
    i_min = i;
  end
  if( vx>v_max )
    v_max = vx;
    i_max = i;
  end
end

out.pass = i_min == 90 & ...
           abs((v_min + 0.5556)/0.5556) < 1e-2 & ...
           i_max == 100 & ...
           abs((v_max - 0.2744)/0.2744) < 1e-2;

if( nargout==0 )
  clear fea out
end
