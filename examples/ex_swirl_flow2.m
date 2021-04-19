function [ fea, out ] = ex_swirl_flow2( varargin )
%EX_SWIRL_FLOW2 2D Axisymmetric swirl flow in step domain.
%
%   [ FEA, OUT ] = EX_SWIRL_FLOW2( VARARGIN ) Axisymmetric swirl for in tubular step
%   region where the inner cylindrical wall is rotating.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       omega       scalar {100}           Angular rotational velocity (of inner wall)
%       sf_u        string {sflag1}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'omega',    100;
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

fea.geom.objects = { gobj_rectangle(0.5,1.5,0,3) gobj_rectangle(1.0,1.5,1.5,3,'R2') };
fea = geom_apply_formula( fea, 'R1-R2' );

fea.grid = gridgen( fea, 'hmax', 0.1, 'fid', fid );


% Equation definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@swirlflow);
  fea.phys.sw.eqn.coef{1,end} = { 1 };
  fea.phys.sw.eqn.coef{2,end} = { 1 };
  fea.phys.sw.sfun            = [ repmat( {opt.sf_u}, 1, 3 ) {opt.sf_p} ];
  fea.phys.sw.bdr.sel = [1 1 5 2 1 1];
  fea.phys.sw.bdr.coef{2,end}{2,4} = opt.omega;
  fea.phys.sw.prop.artstab.ps = isequal(opt.sf_u,opt.sf_p);

  fea = parsephys(fea);
  if( isfield(fea,'constr') )
    fea = rmfield(fea,'constr');
  end

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

  fea.coef = { 'rho', 1 ;
               'miu', 1 ;
               'Fr',  0 ;
               'Fth', 0 ;
               'Fz',  0 };


  % Boundary conditions.
  fea.bdr.d = { 0  0  [] 0  0  0 ;
                0  0  [] opt.omega 0 0 ;
                0  0  0  0  0  0 ;
                [] [] [] [] [] [] };
  fea.bdr.n = cell(size(fea.bdr.d));


end

% Fix pressure at p([r,z]=[ro,h/2]) = 0.
[~,ix_p] = min( sqrt( (fea.grid.p(1,:)-1.5).^2 + (fea.grid.p(2,:)-1.5).^2) );
fea.pnt = struct( 'type',  'constr', ...
                  'index', ix_p, ...
                  'dvar',  'p', ...
                  'expr',  '0' );


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'maxnit', 50, 'fid', fid );


% Postprocessing.
if( opt.iplot )
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', 'isoexpr', 'v' )
end


% Error checking.
out.ref  = [ -6.1 10.5 73 1.25 ];
if( ~got.tol )
  if( opt.sf_u(end) == '2' )
    opt.tol = 0.05;
  else
    opt.tol = 0.3;
  end
end
[u_min,u_max] = minmaxsubd( 'u', fea );
out.val  = [ u_min u_max intsubd('v',fea) intsubd('w',fea) ];
out.pass = mean(abs(out.val-out.ref)./abs(out.ref)) < opt.tol;

if( nargout==0 )
  clear fea out
end
