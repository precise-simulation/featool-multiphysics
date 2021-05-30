function [ fea, out ] = ex_piezoelectric1( varargin )
%EX_PIEZOELECTRIC1 Piezoelectric bending of a beam example.
%
%   [ FEA, OUT ] = EX_PIEZOELECTRIC1( VARARGIN ) Example for piezoelectric bending of a beam.
%
%   References:
%
%   [1] V. Pierfort, Finite Element Modeling of Piezoelectric Active
%   Structures, ULB, Faculty of Applied Sciences, 2000.
%
%   [2] W-S. Hwang, H.C. Park, Finite Element Modeling of Piezoelectric
%   Sensors and Actuators, AIAA Journal, Vol. 31, No. 5, May 1993.
%
%   [3] C-I. Tseng, Electromechanical Dynamics of a Coupled
%   Piezo-electric/Mechanical System Applied to Vibration Control and
%   Distributed Sensing, Univ. of Kentucky, Lexington, July 1989.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       l           scalar {0.12}          Length of beam
%       h           scalar {2e-3}          Height (thickness) of beam
%       delV        scalar {200}           Potential difference
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {20}            Cell resolution
%       sfund       string {sflag1}        Shape function for displacements
%       sfunp       string {sflag1}        Shape function for potential
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'l',        0.12; ...
  'h',        2e-3; ...
  'delV',     100; ...
  'igrid',    0; ...
  'hmax',     20; ...
  'sfund',    'sflag2'; ...
  'sfunp',    'sflag1'; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;
i_impl    = 1;

% Geometry definition.
fea.sdim = { 'x' 'y' };   % Names of space coordinates.
gobj1    = gobj_rectangle( 0, opt.l,  0,       opt.h/2, 'R1' );
gobj2    = gobj_rectangle( 0, opt.l, -opt.h/2, 0,       'R2' );
fea.geom.objects = { gobj1 gobj2 };


% Grid generation.
switch opt.igrid
  case -1
    fea.grid = rectgrid( opt.hmax, opt.hmax, [ 0 opt.l; -opt.h/2 opt.h/2] );
    fea.grid = quad2tri( fea.grid );
  case 0
    fea.grid = rectgrid( opt.hmax, opt.hmax, [ 0 opt.l; -opt.h/2 opt.h/2] );
    fea.grid.s( selcells( fea, 'y<=0' ) ) = 2;
  case 1
    fea.grid = gridgen( fea, 'hmax', opt.h/(opt.hmax/4.6), 'fid', fid );
    fea.grid.s(:) = 1;
end
ind_c = selcells( fea, 'y<=0' );
fea.grid.s( ind_c ) = 2;
ix = ismember( fea.grid.b(1,:), ind_c );
fea.grid.b(3,ix) = fea.grid.b(3,ix) + 4;
[~,~,ib] = unique(fea.grid.b(3,:));
fea.grid.b(3,:) = ib;


% Equation coefficients.
Emod =  2e9;       % Modulus of elasticity
nu   =  0.29;      % Poissons ratio
Gmod =  0.775e9;   % Shear modulus
d31  =  0.22e-10;  % Piezoelectric sn coefficient
d33  = -0.3e-10;   % Piezoelectric sn coefficient
prel = 12;         % Relative electrical permittivity
pvac = 0.885418781762e-11;   % Electrical permittivity of vacuum

% Constitutive relations.
constrel   = [ Emod/(1-nu^2)    nu*Emod/(1-nu^2) 0    ;
               nu*Emod/(1-nu^2) Emod/(1-nu^2)    0    ;
               0                0                Gmod ];
piezoel_st = [ 0 d31 ;
               0 d33 ;
               0 0   ];
piezoel_sn = constrel*piezoel_st;
dielmat_st = [ prel 0    ;
               0    prel ]*pvac ;
dielmat_sn = [ dielmat_st - piezoel_st'*piezoel_sn ];

% Populate coefficient matrices (negative sign due to fem partial integration).
c{1,1} = { constrel(1,1)   constrel(1,3) ;
           constrel(1,3)   constrel(3,3) };
c{1,2} = { constrel(1,3)   constrel(1,2) ;
           constrel(3,3)   constrel(2,3) };
c{2,1} = c{1,2}';
c{2,2} = { constrel(3,3)   constrel(1,3) ;
           constrel(1,3)   constrel(2,2) };
c{1,3} = { piezoel_sn(1,1) piezoel_sn(1,2) ;
           piezoel_sn(3,1) piezoel_sn(3,2) };
if( i_impl )
  for i=1:4
    c{1,3}{i} = [num2str(c{1,3}{i}),'*(2*(y<0)-1)'];
  end
else
  for i=1:4
    fea.expr{i,1} = ['piezoel_sn1',num2str(i)];
    fea.expr{i,2} = { -c{1,3}{i} c{1,3}{i} };
    c{1,3}{i} = ['piezoel_sn1',num2str(i)];
  end
end
c{3,1} = c{1,3}';
c{2,3} = { piezoel_sn(3,1) piezoel_sn(3,2);
           piezoel_sn(2,1) piezoel_sn(2,2) };
if( i_impl )
  for i=1:4
    c{2,3}{i} = [num2str(c{2,3}{i}),'*(2*(y<0)-1)'];
  end
else
  for i=1:4
    fea.expr{i+4,1} = ['piezoel_sn2',num2str(i)];
    fea.expr{i+4,2} = { -c{2,3}{i} c{2,3}{i} };
    c{2,3}{i} = ['piezoel_sn2',num2str(i)];
  end
end
c{3,2} = c{2,3}';
c{3,3} = { dielmat_sn(1,1) dielmat_sn(2,1) ;
           dielmat_sn(2,1) dielmat_sn(2,2) };


% Dependent variable names.
fea.dvar  = { 'u' 'v' 'V' };
n_dvar    = length(fea.dvar);

% Finite element shape functions.
fea.sfun  = { opt.sfund opt.sfund opt.sfunp };

% Define equations.
bilinear_form = [ 2 2 3 3 ;
                  2 3 2 3 ];
for i=1:n_dvar
  for j=1:n_dvar
    fea.eqn.a.form{i,j} = bilinear_form;
    fea.eqn.a.coef{i,j} = c{i,j}(:)';
  end
end

% Source terms (set to zero).
fea.eqn.f.form = { 1 1 1 };
fea.eqn.f.coef = { 0 0 0 };


% Boundary conditions.
n_bdr          = max(fea.grid.b(3,:));   % Number of boundaries.
fea.bdr.d      = cell(n_dvar,n_bdr);
fea.bdr.n      = cell(n_dvar,n_bdr);
i_top          = findbdr( fea, ['y>=',num2str(opt.h/2-sqrt(eps))] );
i_bottom       = findbdr( fea, ['y<=',num2str(-opt.h/2+sqrt(eps))] );
i_left         = findbdr( fea, ['x<=',num2str(sqrt(eps))] );
[fea.bdr.d{3,i_top}]    = deal( opt.delV );   % Set potential to dV on top boundary.
[fea.bdr.d{3,i_bottom}] = deal( 0 );          % Set potential to 0V on bottom boundary.
[fea.bdr.d{1:2,i_left}] = deal( 0 );          % Set displacements to 0 on left boundary.


% Check and parse problem struct.
fea = parseprob( fea );

% Call to stationary solver.
fea.sol.u = solvestat( fea, 'fid', fid );


% Postprocessing.
if( opt.iplot )
  YSCALE = 3;
  axlim  = [0, opt.l, -YSCALE*opt.h/2, YSCALE*opt.h/2];
  DSCALE = 20;

  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'u', 'isoexpr', 'u', 'setaxes', 'off' )
  axis(axlim);
  title('x-displacement')

  subplot(2,2,2)
  postplot( fea, 'surfexpr', 'v', 'isoexpr', 'v', 'setaxes', 'off' )
  axis(axlim);
  title('y-displacement')

  subplot(2,2,3)
  postplot( fea, 'surfexpr', 'V', 'isoexpr', 'V', 'setaxes', 'off' )
  axis(axlim);
  title('Electric potential')

  subplot(2,2,4)
  plotsubd( fea, 'labels', 'off', 'setaxes', 'off' )
  axis(axlim);
  title(['Displacement plot (at ',num2str(DSCALE),' times scale)'])

  ind1 = sub2ind( size(fea.grid.c), fea.grid.b(2,:), fea.grid.b(1,:) );
  p1b  = fea.grid.p( :, fea.grid.c(ind1) );
  ind2 = sub2ind( size(fea.grid.c), mod(fea.grid.b(2,:),size(fea.grid.c,1))+1, fea.grid.b(1,:) );
  p2b  = fea.grid.p( :, fea.grid.c(ind2) );
  up1  = DSCALE*evalexpr( 'u', p1b, fea );
  vp1  = DSCALE*evalexpr( 'v', p1b, fea );
  up2  = DSCALE*evalexpr( 'u', p2b, fea );
  vp2  = DSCALE*evalexpr( 'v', p2b, fea );
  hold on
  for i=1:size(p1b,2)
    plot( [p1b(1,i)+up1(i) p2b(1,i)+up2(i)], [p1b(2,i)+vp1(i) p2b(2,i)+vp2(i)], '-b', 'linewidth', 2 )
  end
  set( gca, 'ytick', [] )
end


% Error checking.
ind_dof_v  = fea.eqn.dofm{2}(:) + fea.eqn.ndof(1);
out.v_max  = max(abs( fea.sol.u(ind_dof_v) ));
v_max_ref  = abs(-3/2*d31*opt.delV*(opt.l/opt.h)^2);
out.err    = abs(v_max_ref - out.v_max)/v_max_ref;
out.pass   = out.err <= 0.1;

if ( nargout==0 )
  clear fea out
end
