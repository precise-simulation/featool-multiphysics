function [ fea, out ] = ex_swirl_flow1( varargin )
%EX_SWIRL_FLOW1 2D Axisymmetric laminar swirl flow.
%
%   [ FEA, OUT ] = EX_SWIRL_FLOW1( VARARGIN ) Axisymmetric swirl for in tubular region
%   where the inner cylindrical wall is rotating. Comparison with analytical solution.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {2}             Density
%       miu         scalar {3}             Molecular/dynamic viscosity
%       omega       scalar {5}             Angular rotational frequency (of inner wall)
%       ri          scalar {0.5}           Inner radius
%       ro          scalar {1.5}           Outer radius
%       h           scalar {3}             Height of cylinder
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'rho',      2;
            'miu',      3;
            'omega',    5;
            'ri',       0.5
            'ro',       1.5;
            'h',        3;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iphys',    1;
            'iplot',    1;
            'tol',      [];
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = {'r' 'z'};

ri = opt.ri;   % Inner radius.
ro = opt.ro;   % Outer radius.
h  = opt.h;    % Height of cylinder.
fea.geom.objects = { gobj_rectangle(ri,ro,-h/2,h/2) };

fea.grid = gridgen( fea, 'hmax', (ro-ri)/8, 'fid', fid );


% Equation definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@swirlflow);
  fea.phys.sw.eqn.coef{1,end} = { opt.rho };
  fea.phys.sw.eqn.coef{2,end} = { opt.miu };
  fea.phys.sw.sfun            = [ repmat( {opt.sf_u}, 1, 3 ) {opt.sf_p} ];
  fea.phys.sw.bdr.sel = [5 1 5 2];
  fea.phys.sw.bdr.coef{2,end}{2,4} = opt.omega*ri;

  fea = parsephys(fea);

else

  opt.sf_u = 'sflag2';
  opt.sf_p = 'sflag1';

  fea.dvar = { 'u', 'v', 'w', 'p' };
  fea.sfun = [ repmat( {opt.sf_u}, 1, 3 ) {opt.sf_p} ];
  c_eqn    = { 'r*rho*u'' - r*miu*(2*ur_r + uz_z  +   wr_z) + r*rho*(u*ur_t + w*uz_t) + r*p_r     = r*Fr - 2*miu/r*u_t + p_t + rho*v*v_t';
               'r*rho*v'' - r*miu*(  vr_r + vz_z) + miu*v_r + r*rho*(u*vr_t + w*vz_t) + rho*u*v_t = r*Fth + miu*(v_r - 1/r*v_t)';
               'r*rho*w'' - r*miu*(  wr_r + uz_r  + 2*wz_z) + r*rho*(u*wr_t + w*wz_t) + r*p_z     = r*Fz';
               'r*ur_t + r*wz_t + u_t = 0' };
  fea.eqn = parseeqn( c_eqn, fea.dvar, fea.sdim );

  fea.coef = { 'rho', opt.rho ;
               'miu', opt.miu ;
               'Fr',  0 ;
               'Fth', 0 ;
               'Fz',  0 };

  % Boundary conditions.
  fea.bdr.d = { []  0 [] 0 ;
                []  0 [] opt.omega*ri ;
                0   0  0 0 ;
                [] [] [] [] };
  fea.bdr.n = cell(size(fea.bdr.d));

  % Fix pressure at p([r,z]=[ro,h/2]) = 0.
  [~,ix_p] = min( sqrt( (fea.grid.p(1,:)-ro).^2 + (fea.grid.p(2,:)-h/2).^2) );
  fea.pnt = struct( 'type',  'constr', ...
                    'index', ix_p, ...
                    'dvar',  'p', ...
                    'expr',  '0' );
end


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', fid );


% Exact (analytical) solution.
a = - opt.omega*ri^2 / (ro^2-ri^2);
b =   opt.omega*ri^2*ro^2 / (ro^2-ri^2);
v_th_ex = @(r,a,b) a.*r + b./r;


% Postprocessing.
if( opt.iplot )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', 'isoexpr', 'v' )

  subplot(1,2,2)
  hold on
  grid on
  r = linspace( ri, ro, 100 );
  v_th = evalexpr( 'v', [r;zeros(1,length(r))], fea );
  plot( r, v_th, 'b--' )
  r = linspace( ri, ro, 10 );
  plot( r, v_th_ex(r,a,b), 'r.' )
  legend( 'Computed solution', 'Exact solution')
  xlabel( 'Radius, r')
  ylabel( 'Angular velocity, v')
end


% Error checking.
if( ~got.tol )
  if( opt.sf_u(end) == '2' )
    opt.tol = 0.01;
  else
    opt.tol = 0.16;
  end
end
r = linspace( ri, ro, 100 );
v_th = evalexpr( 'v', [r;zeros(1,length(r))], fea )';
out.err  = norm( v_th - v_th_ex(r,a,b) );
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
