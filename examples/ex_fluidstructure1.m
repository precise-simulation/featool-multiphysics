function [ fea, out ] = ex_fluidstructure1( varargin )
%EX_FLUIDSTRUCTURE1 Fluid-structure interaction static elastic beam.
%
%   [ FEA, OUT ] = EX_FLUIDSTRUCTURE1( VARARGIN ) Example for fluid-
%   structure interaction flow around a static elastic beam at Re = 20.
%
%   Reference (FSI1):
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
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iplot',    1;
            'tol',      [0.1 0.1 0.1];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


rho   = 1e3;
miu   = 1;
umean = 0.2;
diam  = 0.1;


% Geometry.
fea.sdim = { 'x', 'y' };
gobj1 = gobj_rectangle( 0, 2.5, 0, 0.41, 'R1' );
gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
gobj3 = gobj_rectangle( [0.2], [0.6], [0.2-0.01], [0.2+0.01], 'R2' );
fea.geom.objects = { gobj1 gobj2 gobj3 };
fea.geom = copy_geometry_object( 'C1', fea.geom );
fea.geom = copy_geometry_object( 'R2', fea.geom );
fea.geom = geom_apply_formula( fea.geom, 'R1-C1-R2' );
fea.geom = geom_apply_formula( fea.geom, 'R3-C2' );


% Grid generation.
hmaxb = [ 0.025*ones(1,4) 0.01*ones(1,4) 0.005*ones(1,5) ];
fea.grid = gridgen( fea, 'hmaxb', hmaxb, 'gridgen', 'gridgen2d', 'fid', opt.fid );

p_A = [0.6;0.2];
[~,i_A] = min(sum((fea.grid.p'-repmat(p_A',size(fea.grid.p,2),1)).^2,2));
fea.grid.p(:,i_A) = p_A;

% Equation settings.
fea = addphys( fea, @fluidstructure );
fea.phys.fsi.sfun = { opt.sf_u, opt.sf_u, opt.sf_p, opt.sf_u, opt.sf_u };
fea.phys.fsi.eqn.coef{1,end} = { rho, 1e3 };   % Density.
fea.phys.fsi.eqn.coef{2,end} = { miu,   0 };   % Viscosity.
fea.phys.fsi.eqn.coef{3,end} = { 0,   0.4 };   % Poisson's ratio.
fea.phys.fsi.eqn.coef{4,end} = { 0, 1.4e6 };   % Modulus of elasticity.
fea.phys.fsi.prop.active = [ 1, 0; 1, 0; 1, 0; 0, 1; 0, 1 ];


% Boundary settings.
fea.phys.fsi.bdr.sel = [ 1 3 1 2 1 1 1 1 6 6 -2 -2 -2 ];
fea.phys.fsi.bdr.coef{2,end}{1,4} = ...
[num2str(1.5*umean*4/0.1608),'*y*(0.41-y)*(0.5*(1-cos(pi/2*t))*(t<2)+(t>=2))'];  % Inflow velocity.


% Solver.
fea = parsephys(fea);
fea = parseprob(fea);

[fea.sol.u,fea.sol.t,fea.sol.grid.p] = ...
fsisolve( fea, 'tstep', 0.2, 'tmax', 8, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot(fea,'surfexpr','p')
  subplot(1,2,2)
  postplot(fea,'surfexpr','sqrt(u^2+v^2)')
end


p_Ai = fea.sol.grid.p(:,i_A,end);
ux_A = evalexpr( 'dx', p_Ai, fea );
uy_A = evalexpr( 'dy', p_Ai, fea );


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

  fprintf(fid,'Displacements(A), ux/y = %6g, %6g (Ref: %6g, %6g)\n',ux_A,uy_A,0.0227e-3,0.8209e-3)
  fprintf(fid,'Drag force,       Fd   = %6f (l), %6f (v) (Ref: 14.929)\n',F_d1,F_d2)
  fprintf(fid,'Lift force,       Fl   = %6f (l), %6f (v) (Ref: 1.11905)\n',F_l1,F_l2)
end


% Error checking.
out.uA   = [ux_A,uy_A];
out.Fd   = [F_d1 F_d2];
out.Fl   = [F_l1 F_l2];
out.err  = [abs(out.uA(1)-0.0227e-3)/0.0227e-3, abs(out.uA(2)-0.8209e-3)/0.8209e-3;
            abs(out.Fd-14.29426)/14.29426;
            abs(out.Fl-0.763746)/0.763746];
out.pass = (all(out.err(1,:)<opt.tol(1)))&(out.err(2,2)<opt.tol(2))&(out.err(3,2)<opt.tol(3));
if ( nargout==0 )
  clear fea out
end
