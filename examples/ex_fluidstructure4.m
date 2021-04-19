function [ fea, out ] = ex_fluidstructure4( varargin )
%EX_FLUIDSTRUCTURE4 Fluid-structure interaction - elastic beam.
%
%   [ FEA, OUT ] = EX_FLUIDSTRUCTURE1( VARARGIN ) 3D example for
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
            'igrid',    1;
            'tstep',    0.01;
            'tmax',     1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry.
fea.sdim = { 'x', 'y', 'z' };
gobj1 = gobj_block( 0, 19.5, 0, 12, 0, 0.2, 'B1' );
gobj2 = gobj_block( 3.5, 13.5, 3.5, 8.5, 0, 0.2, 'B2' );
gobj3 = gobj_block( 4.5, 5.5, 5.5, 6.5, 0, 0.2, 'B3' );
gobj4 = gobj_block( 5.5, 9.5, 6-0.03, 6+0.03, 0, 0.2, 'B4' );
fea.geom.objects = { gobj1, gobj2, gobj3, gobj4 };
fea.geom = copy_geometry_object( 'B2', fea.geom );
fea.geom = geom_apply_formula( fea.geom, 'B1-B5' );
fea.geom = copy_geometry_object( 'B4', fea.geom );
fea.geom = geom_apply_formula( fea.geom, 'B2-B6-B3' );


% Grid generation.
if( opt.igrid )
  hmax = [ 0.06 0.8 0.4 ]*0.95;
  fea.grid = gridgen( fea, 'hmax', hmax, 'gridgen', 'gmsh', 'fid', opt.fid );
else
  fea.grid = impexp_gmsh( 'ex_fluidstructure4.msh', 'import' );
end


% Equation and boundary settings.
fea = addphys( fea, @fluidstructure );

rho_f = 1.18e-3;
miu_f = 1.82e-4;
rho_s = 0.1;
nu_s  = 0.25;
E_s   = 2.5e6;

dtol  = num2str(sqrt(eps)*1e3);
ind_s = unique(fea.grid.s(selcells(fea.grid,['(y>=6-0.03-',dtol,').*(y<=6+0.03+',dtol,').*(x>=5.5-',dtol,').*(x<=9.5+',dtol,')'])));
ind_f = setdiff( unique(fea.grid.s), ind_s );

[fea.phys.fsi.eqn.coef{1,end}{ind_f}] = deal(rho_f);
[fea.phys.fsi.eqn.coef{1,end}{ind_s}] = deal(rho_s);
[fea.phys.fsi.eqn.coef{2,end}{ind_f}] = deal(miu_f);
[fea.phys.fsi.eqn.coef{3,end}{ind_s}] = deal(nu_s);
[fea.phys.fsi.eqn.coef{4,end}{ind_s}] = deal(E_s);
fea.phys.fsi.prop.active = zeros(7,length(ind_s)+length(ind_f));
fea.phys.fsi.prop.active(1:4,ind_f) = 1;
fea.phys.fsi.prop.active(5:7,ind_s) = 1;

n_bdr  = max(fea.grid.b(3,:));
i_in   = findbdr( fea, ['x<=',dtol] );
i_out  = findbdr( fea, ['x>=19.5-',dtol] );
i_slip = findbdr( fea, ['y>=12-',dtol,'|y<=',dtol] );
i_per  = [ findbdr( fea, ['z>=0.2-',dtol] ), findbdr( fea, ['z<=',dtol] ) ];
i_fix  = findbdr( fea, ['y>=6-0.03-',dtol,'&y<=6+0.03+',dtol,'&x<=5.5+',dtol] );
i_beam = setdiff( findbdr( fea, ['y>=6-0.03-',dtol,'&y<=6+0.03+',dtol] ), [i_fix,i_per] );
i_ns   = setdiff( findbdr( fea, ['x>=4.5-',dtol,'&x<=5.5+',dtol] ), i_fix );
i_cont = setdiff( 1:max(n_bdr), [i_in,i_out,i_slip,i_per,i_beam,i_ns,i_fix] );

sel = ones(1,n_bdr);
sel( i_in )   =  2;
sel( i_out )  =  3;
sel( i_slip ) =  5;
sel( i_per )  =  5;
sel( i_fix )  =  6;
sel( i_beam ) = -2;
sel( i_cont ) = -1;

fea.phys.fsi.bdr.sel = sel;
fea.phys.fsi.bdr.coef{2,end}{1,i_in} = 51.3;


% Solver.
fea = parsephys(fea);
fea = parseprob(fea);

[fea.sol.u,fea.sol.t,fea.sol.grid.p] = fsisolve( fea, 'tstep', opt.tstep, 'tmax', opt.tmax, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot )
  postplot(fea,'slicex',[],'slicey',[],'sliceexpr','sqrt(u^2+v^2+w^2)')
  rotate3d on
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

out.err  = [ abs((v_min + 0.5)/0.5), abs((v_max - 0.45)/0.45) ];
out.pass = all( out.err < 0.5 );


if( nargout==0 )
  clear fea out
end
