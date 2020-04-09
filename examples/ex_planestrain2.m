function [ fea, out ] = ex_planestrain2( varargin )
%EX_PLANESTRAIN2 Plane strain analysis of a pressure vessel.
%
%   [ FEA, OUT ] = EX_PLANESTRAIN2( VARARGIN ) Benchmark example for plane strain
%   approximation of a pressure vessel (annular cross section with symmetry).
%
%   Reference:
%
%   [1] Barber JR. Solid Mechanics and Its Applications, Intermediate
%   Mechanics of Materials, pp. 461, vol. 175, 2nd ed. Springer, 2011.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.02}          Grid size
%       sfun        string {sflag2}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2020 Precise Simulation, Ltd.


cOptDef = { ...
  'hmax',     0.02;
  'sfun',     'sflag2';
  'use_temp', true;
  'iplot',    1;
  'igrid',    1;
  'tol',      0.1;
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );


E   = 2.1e11;
nu  = 0.3;

ri  = 0.1;
ro  = 0.8;
sri = -5e7;
sro = -2e7;

if( opt.use_temp )
  Ti = 500;
  To = 20;
  K  = (Ti-To)/log(ro/ri);
  al = 1.2e-5;
  T  = @(r) K*log(ro./r) + To;
  Texpr = [num2str(K),'*log(',num2str(ro),'/sqrt(x^2+y^2))+',num2str(To)];
  int_rT = @(r) K*(1/2*r.^2.*log(ro./r) + r.^2/4) + r.^2/2*To;
else
  Ti = 0;
  T0 = 0;
  K  = 0;
  al = 0;
  T  = @(r) 0;
  Texpr = 0;
  int_rT = @(r) 0;
end

% Reference solution.
C = @(r) E*al/(1-nu)./r.^2.*int_rT(r);
B = ( sro - sri + C(ro) - C(ri) )./(1./ro.^2 - 1./ri.^2);
A = sri + C(ri) - B./ri.^2;

ur  = @(r) al*(1+nu)./((1-nu)*r).*int_rT(r) + A*(1-2*nu)*(1+nu)*r/E - (1+nu)*B./(E*r);
sr  = @(r) -E*al/(1-nu)./r.^2.*int_rT(r) + A + B./r.^2;
sth = @(r)  E*al/(1-nu)./r.^2.*int_rT(r) - E*al*T(r)/(1-nu) + A - B./r.^2;
sz  = @(r) nu*(sr(r)+sth(r)) - E*al*(T(r)-To);   % Adding To here.

% Geometry and grid.
fea.sdim = { 'x' 'y' };   % Coordinate names.
gobj = gobj_ellipse( [0 0], [ro], [ro], 'E1' );
fea.geom.objects{1} = gobj;
gobj = gobj_ellipse( [0 0], [ri], [ri], 'E2' );
fea.geom.objects{2} = gobj;
fea.geom = geom_apply_formula( fea.geom, 'E1-E2' );
gobj = gobj_rectangle( [0], [1.2*ro], [0], [1.2*ro], 'R1' );
fea.geom.objects{2} = gobj;
fea.geom = geom_apply_formula( fea.geom, 'CS1&R1' );
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );


% Problem definition.
fea = addphys( fea, @planestrain );
fea.phys.psn.eqn.coef{1,end} = { nu };
fea.phys.psn.eqn.coef{2,end} = { E  };
fea.phys.psn.eqn.coef{6,end} = { al };
fea.phys.psn.eqn.coef{7,end} = { Texpr };
fea.phys.psn.sfun            = { opt.sfun opt.sfun };


% Boundary conditions.
fea.phys.psn.bdr.sel = [ 1, 1, 1, 1 ];
fea.phys.psn.bdr.coef = { '', '', '', {}, { 0, 0, 0, 1; 0, 0, 1, 0 }, [], ...
                          { [num2str(sro),'*nx'], [num2str(sri),'*nx'], '0', '0';
                            [num2str(sro),'*ny'], [num2str(sri),'*ny'], '0', '0' } };


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );


% Error checking.
s_sx = fea.phys.psn.eqn.vars{5,end};
s_sy = fea.phys.psn.eqn.vars{6,end};
s_sz = fea.phys.psn.eqn.vars{7,end};

r = linspace(ri,ro,100);
p = [ r; zeros(size(r)) ];

u    = evalexpr( 'u', p, fea ).';
src  = evalexpr( s_sx, p, fea ).';
sthc = evalexpr( s_sy, p, fea ).';
szc  = evalexpr( s_sz, p, fea ).';
uref   = ur(r);
srref  = sr(r);
sthref = sth(r);
szref  = sz(r);
out.err = [ norm(u-uref)/norm(uref), norm(src-srref)/norm(srref), ...
            norm(sthc-sthref)/norm(sthref), norm(szc-szref)/norm(szref) ];
p = p([2,1],:);
v    = evalexpr( 'v', p, fea ).';
src  = evalexpr( s_sy, p, fea ).';
sthc = evalexpr( s_sx, p, fea ).';
szc  = evalexpr( s_sz, p, fea ).';
out.err = [ out.err;
            norm(v-uref)/norm(uref), norm(src-srref)/norm(srref), ...
            norm(sthc-sthref)/norm(sthref), norm(szc-szref)/norm(szref) ];
out.pass = all( out.err(:) <= opt.tol(:) );


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_sx )
  title( 'Stress in the x-direction' )

  figure
  subplot(2,2,1)
  plot( r, uref, 'b-' )
  hold on
  plot( r, u, 'r.' )
  plot( r, v, 'g.' )
  xlabel('r')
  ylabel('Displacement, u_r')
  grid on

  subplot(2,2,2)
  plot( r, srref, 'b-' )
  hold on
  plot( r, src, 'r.' )
  xlabel('r')
  ylabel('Stress, \sigma_r')
  grid on

  subplot(2,2,3)
  plot( r, sthref, 'b-' )
  hold on
  plot( r, sthc, 'r.' )
  xlabel('r')
  ylabel('Stress, \sigma_\theta')
  grid on

  subplot(2,2,4)
  plot( r, szref, 'b-' )
  hold on
  plot( r, szc, 'r.' )
  xlabel('r')
  ylabel('Stress, \sigma_z')
  grid on
end


if( nargout==0 )
  clear fea out
end
