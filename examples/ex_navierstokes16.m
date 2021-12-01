function [ fea, out ] = ex_navierstokes16( varargin )
%EX_NAVIERSTOKES16 2D Example for stationary flow around a cylinder with an attached beam.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES16( VARARGIN ) Stationary flow
%   around a cylinder with an attached solid beam.
%
%   Reference:
%
%   [1] Hron J. A monolithic FEM/multigrid solver for ALE formulation
%       of fluid structure interaction with application in
%       biomechanics. In H.-J. Bungartz and M. Schäfer, editors,
%       Fluid-Structure Interaction: Modelling, Simulation,
%       Optimisation, LLNCSE. Springer, 2006.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.025}         Grid size
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also EX_NAVIERSTOKES3

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
            'hmax',     0.025;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iplot',    1;
            'tol',      [0.1 0.1];
            'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


rho   = 1e3;
miu   = 1;
umean = 0.2;
diam  = 0.1;


% Geometry.
fea.sdim = { 'x', 'y' };
gobj1 = gobj_rectangle( 0, 2.5, 0, 0.41, 'R1' );
gobj2 = gobj_circle( [0.2 0.2], diam/2, 'C1' );
gobj3 = gobj_rectangle( [0.2], [0.6], [0.2-0.01], [0.2+0.01], 'R2' );
fea.geom.objects = { gobj1 gobj2 gobj3 };
fea.geom = copy_geometry_object( 'C1', fea.geom );
fea.geom = copy_geometry_object( 'R2', fea.geom );
fea.geom = geom_apply_formula( fea.geom, 'R1-C1-R2' );
fea.geom = geom_apply_formula( fea.geom, 'R3-C2' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'gridgen', 'gridgen2d', 'fid', opt.fid );


% Equation settings.
fea = addphys( fea, @navierstokes );
fea.phys.ns.sfun = { opt.sf_u, opt.sf_u, opt.sf_p };
fea.phys.ns.eqn.coef{1,end} = { rho, 0 };   % Density.
fea.phys.ns.eqn.coef{2,end} = { miu, 0 };   % Viscosity.
fea.phys.ns.prop.active = [ 1, 0; 1, 0; 1, 0 ];


% Boundary settings.
fea.phys.ns.bdr.sel = [ 1 3 1 2 ones(1,9) ];
fea.phys.ns.bdr.coef{2,end}{1,4} = ['1.5*',num2str(umean),'/(0.41/2)^2*y*(0.41-y)'];


% Solver.
fea = parsephys(fea);
fea = parseprob(fea);

fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
if ( opt.iplot>0 )
  figure
  subplot(2,1,1)
  postplot( fea, 'surfexpr', s_velm )
  title( 'Velocity field' )
  subplot(2,1,2)
  postplot( fea, 'surfexpr', 'p' )
  title( 'Pressure' )
end


% Calculate benchmark quantities (line integration method).
s_tfx = ['nx*p+',num2str(miu),'*(-2*nx*ux-ny*(uy+vx))'];
s_tfy = ['ny*p+',num2str(miu),'*(-nx*(vx+uy)-2*ny*vy)'];
i_int = [5:8,11:13];   % Integration boundaries.
i_cub = 10;
F_d1  = intbdr(s_tfx,fea,i_int,i_cub,size(fea.sol.u,2),1);
F_l1  = intbdr(s_tfy,fea,i_int,i_cub,size(fea.sol.u,2),1);


% Calculate benchmark quantities (volume integration method).

% Create field 'a' with values one on the cylinder and zero everywhere else.
fea.dvar = [ fea.dvar, {'a'}       ];
fea.sfun = [ fea.sfun, fea.sfun(1) ];
fea      = parseprob(fea);
n_dof    = max(fea.eqn.dofm{1}(:));
[~,ind_gdof] = evalexprdof(0,1,fea,i_int);
u_a      = zeros(n_dof,1);
u_a(ind_gdof) = 1;
fea.sol.u= [fea.sol.u;repmat(u_a,1,size(fea.sol.u,2))];
fea.eqn  = struct;
fea.bdr  = struct;
fea      = parseprob(fea);

s_tfx    = ['ax*p+',num2str(miu),'*(-2*ax*ux-ay*(uy+vx))-(u*ux+v*uy)*a'];
s_tfy    = ['ay*p+',num2str(miu),'*(-ax*(vx+uy)-2*ay*vy)-(u*vx+v*vy)*a'];
F_d2     = intsubd(s_tfx,fea,1,find(fea.grid.s==1),3);
F_l2     = intsubd(s_tfy,fea,1,find(fea.grid.s==1),3);


if( ~isempty(fid) )
  fprintf(fid,'\n\nBenchmark quantities:\n\n')

  fprintf(fid,'Drag force,    Fd = %6f (l), %6f (v) (Ref: 14.929)\n',F_d1,F_d2)
  fprintf(fid,'Lift force,    Fl = %6f (l), %6f (v) (Ref: 1.11905)\n',F_l1,F_l2)
end


% Error checking.
out.Fd   = [F_d1 F_d2];
out.Fl   = [F_l1 F_l2];
out.err  = [abs(out.Fd-14.929)/14.929;
            abs(out.Fl-1.11905)/1.11905];
out.pass = (out.err(1,2)<opt.tol(1))&(out.err(2,2)<opt.tol(2));
if ( nargout==0 )
  clear fea out
end
