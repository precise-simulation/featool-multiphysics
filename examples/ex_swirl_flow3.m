function [ fea, out ] = ex_swirl_flow3( varargin )
%EX_SWIRL_FLOW3 2D Axisymmetric Taylor-Couette (swirl) flow.
%
%   [ FEA, OUT ] = EX_SWIRL_FLOW3( VARARGIN ) Axisymmetric Taylor-Couette swirl flow in
%   a tubular region where the inner cylindrical wall is rotating. Time dependent solution
%   with periodic top and bottom boundaries.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1}             Density
%       miu         scalar {1}             Molecular/dynamic viscosity
%       omega       scalar {300}           Maximum angular velocity (of inner wall)
%       tmax        scalar {3}             Maximum time
%       nstep       scalar {300}           Number of time steps
%       ri          scalar {1.0}           Inner radius
%       ro          scalar {1.5}           Outer radius
%       h           scalar {3}             Height of cylinder
%       sf_u        string {sflag1}        Shape fcn for velocity
%       sf_p        string {sflag1}        Shape fcn for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'rho',      1;
            'miu',      1;
            'omega',    300;
            'tmax',     3;
            'nstep',    300;
            'ri',       1.0
            'ro',       1.5;
            'h',        3;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'iphys',    1;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = {'r' 'z'};

ri = opt.ri;   % Inner radius.
ro = opt.ro;   % Outer radius.
h  = opt.h;    % Height of cylinder.
fea.geom.objects = { gobj_rectangle(ri,ro,-h/2,h/2) };

n = 16;
fea.grid = rectgrid( n, 6*n, [ri ro; -h/2 h/2] );


% Equation definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@swirlflow);
  fea.phys.sw.eqn.coef{1,end} = { opt.rho };
  fea.phys.sw.eqn.coef{2,end} = { opt.miu };
  fea.phys.sw.sfun            = [ repmat( {opt.sf_u}, 1, 3 ) {opt.sf_p} ];
  fea.phys.sw.bdr.sel = [3 1 3 2];
  fea.phys.sw.bdr.coef{2,end}{2,4} = sprintf( '%g*t/%g', opt.omega*ri, opt.tmax );

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
                []  0 [] sprintf( '%g*t/%g', opt.omega*ri, opt.tmax ) ;
                []  0 [] 0 ;
                [] [] [] [] };
  fea.bdr.n = cell(size(fea.bdr.d));

end


% Fix pressure at p([r,z]=[ro,0]) = 0.
[~,ix_p] = min( sqrt( (fea.grid.p(1,:)-ro).^2 + (fea.grid.p(2,:)-0).^2) );
fea.pnt = struct( 'type',  'constr', ...
                  'index', ix_p, ...
                  'dvar',  'p', ...
                  'expr',  '0' );


% Parse and solve problem.
fea = parseprob( fea );
if ( opt.iphys==1 )
  fea.bdr.n{1}{1} = @periodic_vel_bc_1_3;
  fea.bdr.n{2}{1} = @periodic_vel_bc_1_3;
  fea.bdr.n{3}{1} = @periodic_vel_bc_1_3;
else
  [fea.bdr.n{1:3,1}] = deal( @periodic_vel_bc_1_3 );   % Set periodic BCs.
end
fea.sol.u = solvetime( fea, 'ischeme', 1, 'tstep', opt.tmax/opt.nstep, 'tmax', opt.tmax, 'tolchg', inf, 'fid', fid );


% Postprocessing.
if( opt.iplot )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+w^2)', 'isoexpr', 'sqrt(u^2+w^2)', ...
                 'arrowexpr', {'u' 'w'}, 'arrowcolor', 'w', 'arrowspacing', [8 48] )
  title('In-plane velocity')

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'v', 'isoexpr', 'v' )
  title('Azimuthal velocity')
end


out = [];
if( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function l_gen_plots( fea, opt )

ri = opt.ri;
ro = opt.ro;
h  = opt.h;

lw = 2;

close all
f = figure;
axis equal


% Geometry contour lines.
n  = 36 + 1;
th = linspace(0,2*pi,n);
x  = [ ri*cos(th) nan ri*cos(th) nan ro*cos(th) nan ro*cos(th) nan 0 0 nan 0 0 nan ro ro nan -ro -ro ];
y  = [ ri*sin(th) nan ri*sin(th) nan ro*sin(th) nan ro*sin(th) nan ro ro nan -ro -ro nan 0 0 nan 0 0 ];
z  = [ h/2*ones(1,n) nan -h/2*ones(1,n) nan h/2*ones(1,n) nan -h/2*ones(1,n) nan -h/2 h/2 nan -h/2 h/2 nan -h/2 h/2 nan -h/2 h/2 ];
hl = plot3( x, y, z, 'linewidth', lw, 'color', 'k' );

% Inner wall.
n_bands = 4;
verts = [ repmat( ri*cos(th(1:end-1)), 1, 2 ); repmat( ri*sin(th(1:end-1)), 1, 2 ); h/2*ones(1,n-1) -h/2*ones(1,n-1) ]';
faces = [ 1:n-1; 2:n-1 1; n-1+3:2*(n-1) n-1+1 n-1+2; n-1+2:2*(n-1) n-1+1 ]';
cols  = repmat( [1 1 1], n-1, 1 );
ix = mod([1:n-1]-1,(n-1)/n_bands)+1 > (n-1)/n_bands/2;
cols(ix,:) = repmat([0.4 0.4 0.6],sum(ix),1);
hw = patch( 'vertices', verts, 'faces', faces, 'facecolor', 'flat', 'facevertexcdata', cols, 'linestyle', 'none' );

% Arrows
xa = linspace( ri, ro, 8+2 );
xa = xa(2:end-1);
ya = linspace( -h/2,h /2, 48+2 );
ya = ya(2:end-1);
[xx,yy] = meshgrid(xa,ya);
xa = xx(:);
ya = yy(:);
pa = [xa';ya'];
za = ya;
ra = linspace( ri, ro, 8+2 );
ra = ra(2:end-1);
fi = linspace( 0,2 *pi, 16 );

for i=1:size(fea.sol.u,2)

  th = 2*180/pi*i*1e-2;
  rotate( hw, [0 0 1], th/10 )

  % In-plane velocity surface plots.
  h1 = postplot( fea, 'surfexpr', '1e3*sqrt(u^2+w^2)', 'colorbar', 0, 'solnum', i );
  rotate( h1(1), [1 0 0], 90 )
  rotate( h1(1), [0 0 1], -45 )

  h2 = postplot( fea, 'surfexpr', '1e3*sqrt(u^2+w^2)', 'colorbar', 0, 'solnum', i );
  rotate( h2(1), [1 0 0], 90 )
  rotate( h2(1), [0 0 1], 180-45 )

  % Arrows.
  u = evalexpr( 'u', pa, fea, i );
  v = evalexpr( 'v', pa, fea, i );
  w = evalexpr( 'w', pa, fea, i );
  hq = [];
  for j=1:length(fi)
    xa = ra*cos(fi(j));
    ya = ra*sin(fi(j));
    xx = meshgrid(xa,1:48);
    xa = xx(:);
    yy = meshgrid(ya,1:48);
    ya = yy(:);
    hq = [hq quiver3( xa, ya, za, u, v, w, 'color', 'r' )];
  end

  axis off
  drawnow
  print( sprintf('img/img_%03i.jpg',i), '-r300', '-djpeg' )
  delete([h1 h2 hq])
  % pause
end
system('ffmpeg -i img_%03d.jpg -c:v libx264 -vf "fps=25, format=yuv420p" out.mp4')

%------------------------------------------------------------------------------%
function [ A, f, prob ] = periodic_vel_bc_1_3( A, f, prob, i_dvar, i_bdr, solve_step )

if( i_bdr~=1 || solve_step~=1 )   % Only process boundaries 1 ( and 3 ) in
  return                          % the step directly before linear solver.
end
j_bdr = 3;

if( isstruct(A) )
  A = sparse( A.rows, A.cols, A.vals, A.size(1), A.size(2) );
end

bdrm = prob.bdr.bdrm{i_dvar};
ndof = prob.eqn.ndof;
dof_offset = sum( prob.eqn.ndof(1:(i_dvar-1)) );


ix_i     = find( bdrm(3,:)==i_bdr );   % Index to dofs on boundary i_bdr.
[~,itmp] = unique(bdrm(4,ix_i));       % Remove duplicate/shared points (optional).
ix_i     = sort( ix_i(itmp) );         % Sort index.
for i=1:numel(ix_i)
  ix = ix_i(i);
  x_i(i) = evalexpr0( 'r', bdrm(6:end,ix), [], bdrm(1,ix), bdrm(2,ix), prob );
end


ix_j     = find( bdrm(3,:)==j_bdr );   % Index to dofs on boundary j_bdr.
[~,itmp] = unique(bdrm(4,ix_j));       % Remove duplicate/shared points (optional).
ix_j     = sort( ix_j(itmp) );         % Sort index.
ix_ij    = zeros( size(ix_i) );        % Index linking i and j dofs.
for j=1:numel(ix_j)
  ix = ix_j(j);
  x_j = evalexpr0( 'r', bdrm(6:end,ix), [], bdrm(1,ix), bdrm(2,ix), prob );

  [~,ix_ij_j] = min( abs( x_i - x_j ) );
  ix_ij(ix_ij_j) = ix;
end
idofs = bdrm(4,ix_i)  + dof_offset;
jdofs = bdrm(4,ix_ij) + dof_offset;

[~,imin] = min(x_i);
[~,imax] = max(x_i);
ind_rem = [imin imax];
idofs(ind_rem) = [];
jdofs(ind_rem) = [];

A(jdofs,:) = A(idofs,:) + A(jdofs,:);
f(jdofs)   = f(idofs)   + f(jdofs);

A(idofs,:) = 0;
A(sub2ind(size(A),idofs,idofs)) =  1;
A(sub2ind(size(A),idofs,jdofs)) = -1;

f(idofs) = 0;
