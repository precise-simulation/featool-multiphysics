function [ fea, out ] = ex_fluidstructure2( varargin )
%EX_FLUIDSTRUCTURE2 Fluid-structure interaction for an elastic beam.
%
%   [ FEA, OUT ] = EX_FLUIDSTRUCTURE2( VARARGIN ) Example for fluid-
%   structure interaction flow around an elastic beam at Re = 100.
%
%   Reference (FSI2):
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
            'tmax',     15;
            'tol',      [0.1 0.1 0.1];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


rho   = 1e3;
miu   = 1;
umean = 1;
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
hmaxb = [ 0.025*ones(1,4) 0.01*ones(1,4) 0.005*ones(1,5) ]*1;
fea.grid = gridgen( fea, 'hmaxb', hmaxb, 'gridgen', 'gridgen2d', 'fid', opt.fid );

p_A = [0.6;0.2];
[~,i_A] = min(sum((fea.grid.p'-repmat(p_A',size(fea.grid.p,2),1)).^2,2));
fea.grid.p(:,i_A) = p_A;

% Equation settings.
fea = addphys( fea, @fluidstructure );
fea.phys.fsi.sfun = { opt.sf_u, opt.sf_u, opt.sf_p, opt.sf_u, opt.sf_u };
fea.phys.fsi.eqn.coef{1,end} = { rho, 1e4 };   % Density.
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
    fsisolve( fea, 'tstep', 0.01, 'tmax', opt.tmax, 'fid', opt.fid );


% Calculate benchmark quantities (line integration method).
s_tfx = ['nx*p+',num2str(miu),'*(-2*nx*ux-ny*(uy+vx))'];
s_tfy = ['ny*p+',num2str(miu),'*(-nx*(vx+uy)-2*ny*vy)'];
i_int = [5:8,11:13];   % Integration boundaries.
i_cub = 10;

i1 = find(fea.sol.t>=fea.sol.t(end)-1);
i2 = length(fea.sol.t);
i_cnt = 0;
for i=i1:i2
  i_cnt = i_cnt + 1;

  p_Ai = fea.sol.grid.p(:,i_A,i);
  ux_A(i_cnt) = evalexpr( 'dx', p_Ai, fea, i );
  uy_A(i_cnt) = evalexpr( 'dy', p_Ai, fea, i );
  F_d(i_cnt)  = intbdr( s_tfx, fea, i_int, i_cub, i, 1 );
  F_l(i_cnt)  = intbdr( s_tfy, fea, i_int, i_cub, i, 1 );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot(fea,'surfexpr','p')

  figure
  subplot(2,2,1)
  plot( fea.sol.t(i1:i2), ux_A )
  xlabel('time')
  ylabel('x-displacement')

  subplot(2,2,2)
  plot( fea.sol.t(i1:i2), uy_A )
  xlabel('time')
  ylabel('y-displacement')

  subplot(2,2,3)
  plot( fea.sol.t(i1:i2), F_l )
  xlabel('time')
  ylabel('lift force')

  subplot(2,2,4)
  plot( fea.sol.t(i1:i2), F_d )
  xlabel('time')
  ylabel('drag force')
end


% Error checking.
out.t    = fea.sol.t(i1:i2);
out.ux_A = ux_A;
out.uy_A = uy_A;
out.F_d  = F_d;
out.F_l  = F_l;
out.vals = [ min(ux_A), max(ux_A) ;
             min(uy_A), max(uy_A) ;
             min(F_d),  max(F_d)  ;
             min(F_l),  max(F_l)  ];
out.ref  = [ -14.58e-3-12.44e-3, -14.58e-3+12.44e-3 ;
               1.23e-3-80.6e-3,    1.23e-3+80.6e-3 ;
             208.83-73.75,       208.83+73.75 ;
               0.88-234.2,         0.88+234.2 ];
out.err  = abs(out.vals-out.ref)./abs(out.ref);
out.pass = all(out.err(:) < 0.1 | (out.err(:)>=0.1 & out.err(:)<0.5));
if( nargout==0 )
  clear fea out
end
