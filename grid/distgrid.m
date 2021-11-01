function [ p, t, stat ] = distgrid( fcn_dist, fcn_scal, h_init, b_box, p_fix, e_fix, it_max, fid, fit )
%DISTGRID 2D/3D Grid generation using distance functions.
%
%   [ P, T, STAT ] = DISTGRID( FCN_DIST, FCN_SCAL, H_INIT, B_BOX
%                              P_FIX, E_FIX, IT_MAX, FID, FIT )
%
%   DISTGRID is a simple surface (in 2D) and volume (in 3D)
%   grid generation algorithm using distance functions to define
%   geometries.
%
%   FCN_DIST is a function handle to the geometry description that
%   should take evaluation coordinates and points as input. For
%   example fcn_dist = @(p) sqrt(sum(p.^2,2)) - 1; specifies the
%   distance function for a unit circle (both function handles, string
%   function names, and anonymous functions are supported). Similar to
%   FCN_DIST, FCN_SCAL a function describing the desired relative grid
%   size distribution. For example fcn_scal = @(p) ones(size(p,1),1);
%   specifies a uniform distribution where FCN_SCAL evaluates to 1 at
%   all points. H_INIT is a numeric scalar specifying the initial edge
%   lengths, and B_BOX is a 2 by 2 in 2D (or 2 by 3 in 3D) bounding
%   box of the domain (enclosing the zero contour/level set of
%   FCN_DIST). P_FIX optionally specifies a number of points that
%   should always be present (fixed) in the resulting grid. E_FIX can
%   be sets of edge vertex indices to constrain, or alternatively a
%   cell array with function handle to call. IT_MAX sets the maximum
%   number of grid generation iterations allowed (default
%   1200). Finally, FID specifies a file identifies for output
%   (default 1 = terminal output), FIT is an optional function to call
%   every iteration to check for early termination.
%
%   The grid generation function returns the grid point vertices in P,
%   triangulated simplices in T, as well as an optional statistics
%   struct STAT including timings and convergence information.
%
%
%   Input:
%
%      FCN_DIST:  Distance function d(x,y,(z))
%      FCN_SCAL:  Scaled edge length function h(x,y,(z))
%      HINIT:     Initial edge length
%      B_BOX:     Bounding box [xmin,ymin,(zmin); xmax,ymax,(zmax)]
%      P_FIX:     Fixed node positions (n_p_fix x 2/3)
%      E_FIX:     Constrained edges (n_e_fix x 2)
%      IT_MAX:    Maximum number of iterations
%      FID:       Output file id number (default 1 = terminal)
%      FIT:       Optional iteration function call (default none)
%
%   Output:
%
%      P:         Grid vertex/node coordinates (n_p x 2/3)
%      T:         Triangle indices (n_t x 3)
%      STAT:      Grid generation statistics (struct)
%
%
%   Example: (Uniform grid on unit circle)
%      fcn_dist = @(p) sqrt(sum(p.^2,2)) - 1;
%      fcn_scal = @(p) ones(size(p,1),1);
%      [p,t] = distgrid( fcn_dist, fcn_scal, 0.2, [-1,-1;1,1] );
%      patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
%
%   Reference:
%
%   P-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
%   SIAM Review, Volume 46 (2), pp. 329-345, 2004.
%
%   See also GRIDGEN, GRIDGEN_DISTGRID

% Copyright 2013-2021 Precise Simulation, Ltd.
if( ~(nargin || nargout) ),help distgrid, return, end

t0 = tic;
if( nargin<9 )
  fit = [];
end
if( nargin<8 )
  fid = 1;
end
if( nargin<7 )
  it_max = 1200;
end
if( nargin<6 )
  e_fix = [];
end
if( nargin<5 )
  p_fix = [];
end

%------------------------------------------------------------------------------%
% Initialization and grid generation parameters.
%------------------------------------------------------------------------------%
IT_MIN  = 25;                 % Minimum number of iterations.
IT_MINC = 40;                 % Minimum number of iter. after which to call constraint function.
IT_PRT  = 50;                 % Output every IT_PRT iterations.

N_recv  = 2;                  % Number of recovery iteration steps to move points outside back to boundary.
N_dcf   = 30;                 % Frequency of density control checks.
n_sdim  = size(b_box,2);
if( n_sdim==2 )
  dp_tol   = -0.001*h_init;   % Abs point rejection tol (p(dist(p)>=dp0_tol) are rejected).
  dtrm_tol = -0.001*h_init;   % Abs dist tol for tri rejection (t(dist(p_tcent)>=dtrm_tol) are rejected).
  rt_tol   =  0.3;            % Rel fraction of h_init to trigger retriangulation.
  F_scale  =  1.2;            % Rel force scaling factor.
  F_dcf    =  2.0;            % Fraction of L to L_target to allow.
  dp_scale =  0.2;            % Rel fraction of computed new distance to move points in update step.
else
  dp_tol   = -0.1*h_init;
  dtrm_tol = -0.1*h_init;
  rt_tol   =  0.1;
  F_scale  =  1.1;
  F_dcf    =  3.5;
  dp_scale =  0.1;
end
dpc_tol = 0.001*h_init;       % Abs tol for grid point movements during convergence check.
gradeps = sqrt(eps)*h_init;   % Gradient computation offset.
%------------------------------------------------------------------------------%


% Starting grid point distribution, p.
p = l_distgrid_init_p( b_box, h_init );


% Remove points outside geometry (distance function > 0).
ix_keep = -dp_tol > l_distgrid_call_function(fcn_dist,p);
p = p(ix_keep,:);
t = [];
stat = [];
if( isempty(p) )
  return
end

d_keep  = l_distgrid_call_function(fcn_scal,p);
ix_keep = (min(d_keep)./d_keep).^n_sdim > rand(size(p,1),1);
p = p(ix_keep, :);
p_fix = l_distgrid_deduplicate( p_fix );
n_p_fix = size(p_fix,1);
if( ~isempty(p_fix) )
  p = [ p_fix; setdiff(p,p_fix,'rows') ];
end
n_p = size( p, 1 );


% Main grid generation loop.
l_distgrid_message( fid, sprintf('Grid generation:\n\n') )
t1 = tic;
if( it_max<=0 )
  t = l_delaunay_triangulation( p, e_fix );
end
t_tri = toc(t1);
it = 0;
p0 = inf;
n_tri = 0;
n_dcs = 0;
do_break = false;
is_converged = false;
while( it<it_max )
  it = it + 1;

  % Retriangulate if grid points have moved more than tolerance.
  delta_p_max = max( sqrt(sum((p-p0).^2,2)) );
  if( rt_tol*h_init<delta_p_max )
    [p,t,ep,n_tri,t_tri] = l_distgrid_triangulate( p, fcn_dist, dtrm_tol, e_fix, n_p_fix, it, IT_MINC, n_tri, t_tri );
    n_p = size(p,1);
  end


  % Determine forces and move grid points.
  p1 = p(ep(:,1),:);
  p2 = p(ep(:,2),:);
  edges = p1 - p2;
  len_edges = sum(edges.*edges,2).^(1/2);
  h_edges = l_distgrid_call_function( fcn_scal, 0.5*( p1 + p2 ) );
  len_target = F_scale*h_edges * ( sum(len_edges.^n_sdim) / sum(h_edges.^n_sdim) );
  len_target = len_target.^(1/n_sdim);

  % Remove points not fulfilling criteria.
  ix_crit = F_dcf*len_edges < len_target;
  if( mod(it,N_dcf)==0 && any(ix_crit) )
    n_dcs = n_dcs + 1;
    ind_keep = ep(~ix_crit,:);
    ind_keep = unique([[1:n_p_fix].';ind_keep(:)]);
    p = p(ind_keep,:);
    n_p = size(p,1);
    p0  = inf;
    continue;
  end

  % Compute grid point movements.
  F_edge = repmat( max( [len_target - len_edges, zeros(size(len_target))], ...
                        [], 2 ) ./ len_edges, 1,n_sdim) .* edges;
  delta_p = [];
  for i=1:n_sdim
    delta_p = [ delta_p, ...
                accumarray(ep(:),[F_edge(:,i); -F_edge(:,i)],[n_p,1]) ];
  end
  delta_p(1:n_p_fix,:) = 0;
  delta_p = dp_scale * delta_p;
  p = delta_p + p;


  % Move grid points with distance>0 back to geometry.
  TOL_GRADNM = eps*1e3;
  for jt=1:N_recv

    dist = l_distgrid_call_function( fcn_dist, p );
    ix = dist > 0;
    ix(1:n_p_fix) = 0;

    if( any(ix) )
      grad_dist   = zeros(sum(ix),n_sdim);
      for i=1:n_sdim
        doff = zeros(1,n_sdim);
        doff(i) = gradeps;

        dist_offset_i = l_distgrid_call_function( fcn_dist, p(ix,:)+ones(sum(ix),1)*doff );
        grad_dist(:,i) = ( dist_offset_i - dist(ix) )/gradeps;
      end
      gradnm = sum( grad_dist.^2, 2 );

      jx = gradnm < TOL_GRADNM;
      if( any(jx) )   % Central difference if one sided differentiation fails.
        jx  = find(jx);
        iix = find(ix);
        for i=1:n_sdim
          doff = zeros(1,n_sdim);
          doff(i) = gradeps;

          dist_offset_i = l_distgrid_call_function( fcn_dist, p(iix(jx),:)+ones(length(jx),1)*doff );
          dist_offset_j = l_distgrid_call_function( fcn_dist, p(iix(jx),:)-ones(length(jx),1)*doff );
          grad_dist(jx,i) = ( dist_offset_i - dist_offset_j )/gradeps*0.5;
        end
        gradnm = sum( grad_dist.^2, 2 );
      end
      gradnm(gradnm<TOL_GRADNM) = realmax;

      p(ix,:) = p(ix,:) - (dist(ix)./gradnm*ones(1,n_sdim) .* grad_dist);
    end

  end


  % Statistics and terminal output.
  delta_p_max = abs( max( [sqrt(sum(delta_p(dist<dp_tol,:).^2,2)); -inf] ) );
  if( it==1 || ~mod(it,IT_PRT) )
    s = sprintf( 'Iteration %4i: %i vertices, %i cells, max(delta_p) = %g\n', it, size(p,1), size(t,1), delta_p_max );
    l_distgrid_message( fid, s )
    if( ~isempty(fit) )
      do_break = l_distgrid_call_function( fit, it );
    end
  end


  % Check for convergence.
  if( (it>IT_MIN && delta_p_max<dpc_tol) || size(t,1)<=2 || it>it_max || do_break )
    if( delta_p_max<dpc_tol )
      is_converged = true;
    end
    break;
  end

end


% Correct and check final grid.
[p,t,td] = l_correct_grid( p, t, fcn_dist, e_fix, dtrm_tol, false, n_sdim==2 );
t_tri = t_tri + td;


% Statistics.
t_tot = toc(t0);
s = sprintf( 'Iteration %4i: %i vertices, %i cells, max(delta_p) = %g\n', it, size(p,1), size(t,1), delta_p_max );
l_distgrid_message( fid, s )
if( nargout>=3 )
  stat.conv = is_converged;
  stat.nit  = it;
  stat.ntri = n_tri;
  stat.ndcs = n_dcs;
  stat.dpmx = max(sqrt(sum(delta_p(dist<-dp_tol,:).^2,2)));
  stat.dpmn = mean(sqrt(sum(delta_p(dist<-dp_tol,:).^2,2)));
  stat.ttot = t_tot;
  stat.ttri = t_tri;
end


%------------------------------------------------------------------------------%
function [ p ] = l_distgrid_init_p( b_box, h_init );

n_sdim = size(b_box,2);
p_init = cell(1,n_sdim);
for i=1:n_sdim
  scal = 1;
  if( n_sdim==2 && i==2 )
    scal = sqrt(3)/2;
  end
  p_init{i} = b_box(1,i):(h_init*scal):b_box(2,i);
end

if( n_sdim==2 )
  [p_tmp1,p_tmp2] = ndgrid( p_init{:} );
  p_tmp1(:,2:2:end) = p_tmp1(:,2:2:end) + h_init/2;
  p = [ p_tmp1(:), p_tmp2(:) ];
else   % n_sdim==3
  [p_tmp1,p_tmp2,p_tmp3] = ndgrid( p_init{:} );
  p = [ p_tmp1(:), p_tmp2(:), p_tmp3(:) ];
end

%------------------------------------------------------------------------------%
function [ p, t, edge_pairs, n_tri, t_tri ] = l_distgrid_triangulate( p, fcn_dist, dtrm_tol, e_fix, n_p_fix, it, IT_MINC, n_tri, t_tri )

n_tri = n_tri + 1;

[p,t,td] = l_distgrid_delaunay( p, fcn_dist, e_fix, dtrm_tol );
t_tri = t_tri + td;

n_sdim = size(p,2);
if( iscell(e_fix) && it>IT_MINC )
  [p,t,do_retri] = l_distgrid_call_function( e_fix, p, t, n_sdim, 1:n_p_fix );
  if( do_retri )
    n_tri = n_tri + 1;
    [p,t,td] = l_distgrid_delaunay( p, fcn_dist, e_fix, dtrm_tol );
    t_tri = t_tri + td;
  end
end
p0  = p;
n_p = size(p,1);

e = [ t(:,[1,2]); t(:,[2,3]); t(:,[3,1]) ];
if( n_sdim==3 )
  e = [ e; t(:,[1,4]); t(:,[2,4]); t(:,[3,4]) ];
end
e = sort(e,2);
e_max = max(e(:));
if( e_max*(e_max+1)<realmax )
  ecomp = (e_max+1)*e(:,1) + e(:,2);
  [tmp,ind] = unique( ecomp );
  edge_pairs = e(ind,:);
else
  edge_pairs = unique( e, 'rows' );
end

%------------------------------------------------------------------------------%
function [ p, t, td ] = l_distgrid_delaunay( p, fd, e_fix, dtrm_tol, check_vol )

if( nargin<5 )
  check_vol = true;
end
if( nargin<3 )
  e_fix = [];
end


AV_TOL = prod(abs(max(p,[],1)-min(p,[],1)))*1e-9;   % Minimum accepted absolute area/volume.

[is_nan,tmp] = find( isnan(p) );
p(is_nan,:)  = [];
p = l_distgrid_deduplicate( p );


% Generate triangulation for grid points p.
t1 = tic;
t = l_delaunay_triangulation( p, e_fix );
td = toc(t1);


% Calculate grid cell centers.
pc = [];
for i=1:size(p,2)
  pc = [ pc, mean(reshape(p(t,i),size(t)),2) ];
end

% Remove cells with center outside region.
dist = l_distgrid_call_function(fd,pc);
t = t(dist<dtrm_tol,:);


% Reorient cells.
av = l_grid_cell_volume( p, t );
ix_flip = av<0;
t(ix_flip,[1,2]) = t(ix_flip,[2,1]);


% Remove cells with volume < AV_TOL.
if( check_vol )
  t(abs(av)<AV_TOL,:) = [];
end


if( isempty(t) )
  t = l_delaunay_triangulation( p, e_fix );
end

%------------------------------------------------------------------------------%
function [ t ] = l_delaunay_triangulation( p, c )

if( nargin<2 )
  c = [];
end

IS_WARN       = false;
USE_DELAUNAYN = false;
IS_CONSTR     = isnumeric(c) & ~isempty(c);

if( size(p,2)==3 && USE_DELAUNAYN )
  t = delaunayn( p );
else
  if( ~isempty(c) && IS_CONSTR && exist('DelaunayTriangulation') && size(p,2)==2 )
    try
      if( ~IS_WARN )
        warning('off','MATLAB:delaunayTriangulation:ConsSplitPtWarnId')
      end
      t = delaunayTriangulation( p, c );
      if( ~IS_WARN )
        warning('on','MATLAB:delaunayTriangulation:ConsSplitPtWarnId')
      end
    catch
      t = delaunay( p );
    end

  elseif( exist('delaunayTriangulation') )

    t = delaunayTriangulation( p );

  elseif( exist('DelaunayTri') )

    if( size(p,2)==3 )
      t = DelaunayTri( p(:,1), p(:,2), p(:,3) );
    else
      t = DelaunayTri( p(:,1), p(:,2) );
    end

  elseif( size(p,2)==3 && exist('delaunay3') && ...
          ~exist('OCTAVE_VERSION','builtin') )

    t = delaunay3( p(:,1), p(:,2), p(:,3) );

  else

    t = delaunay( p );

  end
end

if( isa(t,'DelaunayTri') )
  t = t.Triangulation;
end
if( isa(t,'delaunayTriangulation') )
  t = t.ConnectivityList;
end

%------------------------------------------------------------------------------%
function [ p, t, td, ind_p ] = l_correct_grid( p, t, fd, e_fix, dtrm_tol, check_vol, check_qual )
% Remove duplicated/unused nodes and correct element orientation.

if( nargin<7 )
  check_qual = false;
end
if( nargin<6 )
  check_vol = true;
end
if( nargin>=2 && (isempty(p) || isempty(t)) )
  ind_p = 1:size(p,1);
  return
end


P_TOL = eps*1024;
[p,ix,ind_p_orig] = l_distgrid_deduplicate( p, P_TOL );


if( nargin>=2 )
  t = ind_p_orig(t);

  % Final triangulation.
  [p,t,td] = l_distgrid_delaunay( p, fd, e_fix, dtrm_tol, check_vol );

  % Calculate simplex centers.
  pc = [];
  for i=1:size(p,2)
    pc = [ pc, mean(reshape(p(t,i),size(t)),2) ];
  end

  % Remove simplices with center outside region.
  dist = l_distgrid_call_function(fd,pc);
  t = t(dist<dtrm_tol,:);

  % Remove cells with quality < QUAL_TOL.
  if( check_qual )
    a = sqrt(sum((p(t(:,1),:) - p(t(:,2),:)).^2,2));
    b = sqrt(sum((p(t(:,2),:) - p(t(:,3),:)).^2,2));
    c = sqrt(sum((p(t(:,3),:) - p(t(:,1),:)).^2,2));
    r = 0.5*sqrt((b+c-a).*(a+b-c).*(c+a-b)./(a+b+c));
    R = a.*b.*c./sqrt((a+b+c).*(c+a-b).*(b+c-a).*(a+b-c));
    ix = R<eps;
    R(ix) = 1;
    q = 2*r./R;

    QUAL_TOL = 0.03;
    t(q<QUAL_TOL,:) = [];
  end

  % Remove unused nodes.
  [ind_p,ix1,jx1] = unique( t );
  t = reshape( jx1, size(t) );
  p = p( ind_p, : );
  ind_p = ix( ind_p );

end

%------------------------------------------------------------------------------%
function [ v ] = l_grid_cell_volume( p, t )

d12 = p(t(:,2),:) - p(t(:,1),:);
d13 = p(t(:,3),:) - p(t(:,1),:);
if( size(p,2)==2 )
  v = ( d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1) )/2;
else
  d14 = p(t(:,4),:) - p(t(:,1),:);
  v = dot( cross(d12,d13,2), d14, 2 )/6;
end

%------------------------------------------------------------------------------%
function [ b, i, j ] = l_distgrid_deduplicate( a, atol )

if( isempty(a) )
  b = []; return;
end
if( nargin>=2 )
  s = atol;
else
  TOL = 1e-6;
  s = TOL*max(max(a)-min(a));
end
[c,k] = sortrows(s*round(a/s));
ix = any(c(1:size(c,1)-1,:)~=c(2:size(c,1),:),2);
j(k) = cumsum([1;ix]);
i = k([1;find(ix)+1]);
if( nargout>2 )
  [i,jj] = sort(i);
  kk(jj) = 1:numel(jj);
  j = kk(j).';
else
  i = sort(i);
end
b = a(i,:);

%------------------------------------------------------------------------------%
function [ varargout ] = l_distgrid_call_function( fun, varargin )

if( isa(fun,'function_handle') )

  varargout = cell(1,max(1,nargout(fun)));
  [varargout{:}] = fun( varargin{:} );

elseif( iscell(fun) && (isa(fun{1},'function_handle') || ischar(fun{1})) )
  args = fun(2:end);
  fun  = fun{1};
  if( ischar(fun) )
    fun = str2func(fun);
  end

  empty_pos = find(cellfun(@isempty,args));
  if( ~isempty(empty_pos) )
    for i=1:length(varargin)
      if( i<=length(empty_pos) )
        args{empty_pos(i)} = varargin{i};
      else
        args = [ args, varargin{i} ];
      end
    end
  else
    args = [ varargin, args ];
  end

  varargout = cell(1,max(1,nargout(fun)));
  [varargout{:}] = fun( args{:} );

end

if( nargout>0 && ~iscell(varargout) )
  varargout = { varargout };
end

%------------------------------------------------------------------------------%
function l_distgrid_message( fid, s )

if( isscalar(fid) && isnumeric(fid) && fid>0 )
  if( ~any(double(s(end))==[10,13]) )
    s = [s,char(10)];
  end
  fprintf( fid, s );
elseif( isa(fid,'function_handle') )
  fid( s );
end
