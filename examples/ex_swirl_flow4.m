function [ fea, out ] = ex_swirl_flow4( varargin )
%EX_SWIRL_FLOW4 2D Rotating swirling flow around a disk.
%
%   [ FEA, OUT ] = EX_SWIRL_FLOW4( VARARGIN ) Axisymmetric swirling flow around a
%   rotating disk immersed in a container.
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
            'iphys',    1;
            'iplot',    1;
            'tol',      0.01;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = {'r' 'z'};

fea.geom.objects = { gobj_polygon([0 0;.03 0;.03 .05;.002 .05;.002 .02;.008 .018;.0085 .013;0 .013]) };


hmax = 2e-3;
fea.grid = gridgen( fea, 'hmax', hmax, 'fid', fid );


% Equation definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@swirlflow);
  fea.phys.sw.eqn.coef{1,end} = { 1.2e3 };
  fea.phys.sw.eqn.coef{2,end} = { 2.3e-3 };
  fea.phys.sw.sfun            = [ repmat( {opt.sf_u}, 1, 3 ) {opt.sf_p} ];
  fea.phys.sw.bdr.sel = [1 1 3 2 2 2 2 5];
  [fea.phys.sw.bdr.coef{2,end}{2,4:7}] = deal('r*pi/6');

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

  fea.coef = { 'rho', 1.2e3 ;
               'miu', 2.3e-3 ;
               'Fr',  0 ;
               'Fth',  0 ;
               'Fz',  0 };


  % Boundary conditions.
  n_bdr = max(fea.grid.b(3,:));
  fea.bdr.d = cell(4,n_bdr);
  [fea.bdr.d{1,[1 2 5:7 8]}] = deal(0);
  [fea.bdr.d{3,[1 2 3 4:7]}] = deal(0);
  [fea.bdr.d{2,[1 2]}] = deal(0);
  [fea.bdr.d{2, 4:7 }] = deal('r*pi/6');
  [fea.bdr.d{1:3,3}] = deal([]);
  fea.bdr.n = cell(size(fea.bdr.d));


  % Fix pressure at p([r,z]=[0.03,0.05]) = 0.
  [~,ix_p] = min( sqrt( (fea.grid.p(1,:)-0.03).^2 + (fea.grid.p(2,:)-0.05).^2) );
  fea.pnt = struct( 'type',  'constr', ...
                    'index', ix_p, ...
                    'dvar',  'p', ...
                    'expr',  '0' );
end


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', fid );


% Postprocessing.
if( opt.iplot )
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', 'isoexpr', 'v', 'isolev', 15, 'isocolor', 'w' )
end


% Error checking.
[~,U_max] = minmaxsubd( 'sqrt(u^2+v^2+w^2)', fea );
out.err  = norm( U_max - 4.445e-3 )/4.445e-3;
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
