function [ fea, out ] = ex_navierstokes14( varargin )
%EX_NAVIERSTOKES14 Axisymmetric flow in a pipe due to pressure difference.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES14( VARARGIN ) Sets up and solves
%   stationary Poiseuille flow in an axisymmetric pipe driven by a
%   pressure difference. The flow profile is constant and should
%   assume a parabolic profile u(r)=dp/dx/2/miu*(r+h/2)*(h/2-r).
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {0.1}           Density
%       miu         scalar {0.2}           Molecular/dynamic viscosity
%       dp          scalar {0.3}           Pressure difference between in and out-flow
%       h           scalar {0.5}           Pipe/channel height/width
%       l           scalar {2.5}           Pipe/channel length
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.04}          Max grid cell size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'rho',      0.1;
            'miu',      0.2;
            'dp',       0.3;
            'h',        0.5;
            'l',        2.5;
            'igrid',    1;
            'hmax',     0.5/10;
            'sf_u',     'sflag2';
            'sf_p',     'sflag1';
            'iphys',    1;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model parameters.
rho       = opt.rho;     % Density.
miu       = opt.miu;     % Molecular/dynamic viscosity.
dp        = opt.dp;      % Pressure differetial.
% Geometry and grid parameters.
h         = opt.h;       % Height/width of pipe.
l         = opt.l;       % Length of pipe.
% Discretization parameters.
sf_u      = opt.sf_u;    % FEM shape function type for velocity.
sf_p      = opt.sf_p;    % FEM shape function type for pressure.


% Geometry definition.
gobj = gobj_rectangle( 0, h/2, 0, l );
fea.geom.objects = { gobj };
fea.sdim = { 'r' 'z' };   % Coordinate names.


% Grid generation.
fea.grid = rectgrid(round(l/opt.hmax),round(h/2/opt.hmax),[0 h/2;0 l]);
if( opt.igrid~=0 )
  fea.grid = quad2tri( fea.grid );
end


% Boundary conditions.
dtol      = opt.hmax/2;
i_inflow  = findbdr( fea, ['z<',num2str(dtol)] );     % Inflow boundary number.
i_outflow = findbdr( fea, ['z>',num2str(l-dtol)] );   % Outflow boundary number.
i_axis    = findbdr( fea, ['r<',num2str(dtol)] );     % Symmetry axis boundary number.


% Problem definition.
fea = addphys(fea,@navierstokes,'ns',true);   % Add Navier-Stokes equations physics mode.
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
fea.phys.ns.sfun = { sf_u sf_u sf_p };   % Set shape functions.


fea.phys.ns.bdr.sel([i_inflow i_outflow]) = 4;
fea.phys.ns.bdr.coef{4,end}{3,i_inflow}   = dp;         % Set inflow pressure.
fea.phys.ns.bdr.coef{4,end}{3,i_outflow}  = 0;          % Set outflow pressure.
fea.phys.ns.bdr.sel(i_axis)               = 5;

fea       = parsephys(fea);                 % Check and parse physics modes.
fea       = parseprob( fea );               % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'fid', fid );   % Call to stationary solver.


% Postprocessing.
s_velm = 'sqrt(u^2+v^2)';
s_refsol = [num2str(dp),'/',num2str(l),'/2/',num2str(miu),'*(r+',num2str(h),'/2)*(',num2str(h),'/2-r)'];   % Definition of velocity profile.
p = [ linspace(0,h/2,25); l/2*ones(1,25) ];
u = evalexpr( s_velm, p, fea );
u_ref = evalexpr( s_refsol, p, fea );
if ( opt.iplot>0 )
  figure
  subplot(1,3,1)
  postplot(fea,'surfexpr',s_velm,'evaltype','exact')
  title('Velocity field')

  subplot(1,3,2)
  postplot(fea,'surfexpr','p','evaltype','exact')
  title('Pressure')

  subplot(1,3,3)
  plot( u, p(1,:) )
  hold on
  plot( u_ref, p(1,:), 'k.' )
  legend( 'Computed solution', 'Reference solution','Location','West')
  title( 'Velocity profile at z=l/2' )
  ylabel( 'r' )
end


% Error checking.
err = sqrt(sum((u-u_ref).^2)/sum(u_ref.^2));
if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.pass = err<0.05;
if ( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ A, f, prob ] = periodic_pressure_bc_2_4( A, f, prob, i_dvar, i_bdr, solve_step )

if( i_bdr~=2 || solve_step~=1 )   % Only process boundaries 2 ( and 4 ) in
  return                          % the step directly before linear solver.
end
j_bdr = 4;

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
  y_i(i) = evalexpr0( 'y', bdrm(6:end,ix), [], bdrm(1,ix), bdrm(2,ix), prob );
end


ix_j     = find( bdrm(3,:)==j_bdr );   % Index to dofs on boundary j_bdr.
[~,itmp] = unique(bdrm(4,ix_j));       % Remove duplicate/shared points (optional).
ix_j     = sort( ix_j(itmp) );         % Sort index.
ix_ij    = zeros( size(ix_i) );        % Index linking i and j dofs.
for j=1:numel(ix_j)
  ix = ix_j(j);
  y_j = evalexpr0( 'y', bdrm(6:end,ix), [], bdrm(1,ix), bdrm(2,ix), prob );

  [~,ix_ij_j] = min( abs( max(y_i) - y_i - y_j ) );   % Reversed order.
  ix_ij(ix_ij_j) = ix;
end
idofs = bdrm(4,ix_i)  + dof_offset;
jdofs = bdrm(4,ix_ij) + dof_offset;

A(jdofs,:) = A(idofs,:) + A(jdofs,:);
f(jdofs)   = f(idofs)   + f(jdofs);

A(idofs,:) = 0;
A(sub2ind(size(A),idofs,idofs)) =  1;
A(sub2ind(size(A),idofs,jdofs)) = -1;

global dp
f(idofs) = dp;
