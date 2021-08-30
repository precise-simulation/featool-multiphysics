function [ fea, out ] = ex_spanner( varargin )
%EX_SPANNER Gmsh grid import and stress calculation of a spanner.
%
%   [ FEA, OUT ] = EX_SPANNER( VARARGIN ) Example to import a STL CAD
%   geometry or Gmsh grid and calculate displacements on a fixed
%   spanner. The load force may be distributed in the tangential load
%   direction with the force fraction parameter FRAC.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {190e3}         Modulus of elasticity [N/mm^2]
%       nu          scalar {0.29}          Poissons ratio
%       force       scalar {1000}          Load force [N]
%       frac        scalar {0}             Fraction of stress against pulling direction
%       sfun        string {sflag1}        Shape function for displacements
%       igeom       scalar 0/{1}           Import STL geometry (=1)/import Gmsh grid (=0)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'E',        190e3;
            'nu',       0.29;
            'force',    1000;
            'frac',     0;
            'sfun',     'sflag1';
            'igeom',    1;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Define scaling factor m to mm.
USE_METERS = true;
s = double( ~USE_METERS + USE_METERS*1e-3 );


% Import (and scale) grid.
fea.sdim = {'x','y','z'};
if( opt.igeom )
  fea.geom = impexp_stl( 'spanner1.stl', 'import', [], 'solid', 2, 'extend', 0 );
  fea.grid = gridgen( fea, 'hmax', 5, 'gridgen', 'robust', 'fid', opt.fid );
else
  fea.grid = impexp_gmsh( 'spanner.msh', 'import', [], [], opt.fid );
  fea.geom.objects{1} = gobj_grid( fea.grid );
end
fea.grid.p = fea.grid.p*s;


% Add linear elasticity physics mode and define material parameters.
fea = addphys( fea, @linearelasticity );
fea.phys.el.sfun            = { opt.sfun opt.sfun opt.sfun };
fea.phys.el.eqn.coef{1,end} = { opt.nu };
fea.phys.el.eqn.coef{2,end} = { opt.E/s^2 };


% Set all boundaries to no load per default.
n_bdr  = max(fea.grid.b(3,:));
bc_sel = cell(3,n_bdr);
[bc_sel{:}] = deal(0);


% Fix all displacements on mandible boundaries.
if( opt.igeom )
  i_fix   = [13 14];
  i_force = [3 1];
  faxis   = [1 3];
else   % Read mesh.
  i_fix   = [11 14];
  i_force = [1 4];   % Right (pos x)/bottom (neg y) boundary.
  faxis   = [1 3];
end
[bc_sel{:,i_fix}] = deal(1);
fea.phys.el.bdr.coef{5} = bc_sel;


% Apply force for x > 140 mm.
force = opt.force/(6*s*80*s);
fea.phys.el.bdr.coef{7}{faxis(1),i_force(1)} = ['-',num2str((1-opt.frac)*force),'*(y>140*',num2str(s),')'];
fea.phys.el.bdr.coef{7}{faxis(2),i_force(2)} = [num2str(opt.frac*force),'*(y>140*',num2str(s),')'];


% Parse and solve problem.
fea = parsephys(fea);
fea = parseprob(fea);
fea.sol.u = solvestat( fea, 'fid', opt.fid, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')) );


% Postprocessing.
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', ['sqrt(u^2+v^2+w^2)/',num2str(s)] )
  view([30 20])
  title('Total displacement (mm)')

  subplot(1,2,2)
  DSCALE = 5;
  dp = zeros(size(fea.grid.p));
  for i=1:3
    dp(i,:) = DSCALE*evalexpr( fea.dvar{i}, fea.grid.p, fea );
  end
  fea_disp.grid   = fea.grid;
  fea_disp.grid.p = fea_disp.grid.p + dp;
  plotgrid( fea_disp )
  title(['Displacement plot'])
  view([30 20])
end


% Error checking.
u = fea.sol.u(unique(fea.eqn.dofm{1}(:)))/s;
v = fea.sol.u(unique(fea.eqn.ndof(1)+fea.eqn.dofm{2}(:)))/s;
w = fea.sol.u(unique(sum(fea.eqn.ndof(1:2))+fea.eqn.dofm{3}(:)))/s;
out.disp = sqrt( u.^2 + v.^2 + w.^2 );
out.pass = nan;
if( ~(got.frac || got.E || got.nu || got.force) )
  out.pass = abs( max(out.disp) - 2.5 )/2.5 < 0.1;
end


if( nargout==0 )
  clear fea out
else
  fea.grid.p = fea.grid.p/s;
end
