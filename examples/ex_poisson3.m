function [ fea, out ] = ex_poisson3( varargin )
%EX_POISSON3 2D Poisson equation example on a unit square.
%
%   [ FEA, OUT ] = EX_POISSON3( VARARGIN ) Poisson equation on a [0..1]^2 unit square.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Max grid cell size
%       sfun        string {sflag1}        Shape function
%       iphys       scalar 0/{1}           Use physics mode to define problem    (=1)
%                                          or directly define fea.eqn/bdr fields (=0)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'igrid',    0; ...
  'hmax',     0.1; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'icub',     2; ...
  'iplot',    1; ...
  'tol',      0.01;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.geom.objects = { gobj_rectangle() };


% Grid generation.
switch opt.igrid
  case -2
    n  = round(1/opt.hmax);
    fea.grid = rectgrid(n,n,[0 1;0 1]);
    ix = setdiff( 1:(n+1)^2, [1:(n+1) (n+1):(n+1):(n+1)^2 (n+1)^2:-1:(n+1)^2-(n+1)+1 (n+1)^2-(n+1)+1:-(n+1):1 ] );
    h  = 0.3*1/n;
    fea.grid.p(1,ix) = fea.grid.p(1,ix) + h*(2*rand(1,numel(ix)) - 1);
    fea.grid.p(2,ix) = fea.grid.p(2,ix) + h*(2*rand(1,numel(ix)) - 1);
  case -1
    fea.grid = rectgrid(round(1/opt.hmax),round(1/opt.hmax),[0 1;0 1]);
    fea.grid = quad2tri(fea.grid);
  case 0
    fea.grid = rectgrid(round(1/opt.hmax),round(1/opt.hmax),[0 1;0 1]);
  case 1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
end
n_bdr = max(fea.grid.b(3,:));           % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };                % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { 1 };   % Set source term coefficient.
  fea.phys.poi.bdr.coef{1,end} = repmat({0},1,n_bdr);   % Set Dirichlet boundary coefficient to zero.
  fea = parsephys(fea);                 % Check and parse physics modes.

  if( any(strcmp(opt.sfun,{'sf_tri_H3','sf_quad_H3'})) )   % Prescribed derivatives at end points for Hermite elements.
    fea.bdr.d = {{ 0 0 0 0 ;
                   0 '-ex_poisson3_derivative(y)' 0 'ex_poisson3_derivative(y)' ;
                   'ex_poisson3_derivative(x)' 0 '-ex_poisson3_derivative(x)' 0 }};
  end

else

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun };             % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2 3;2 3] };       % First row indicates test function space   (2=x-derivative + 3=y-derivative),
                                        % second row indicates trial function space (2=x-derivative + 3=y-derivative).
  fea.eqn.a.coef = { 1 };               % Coefficient used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { 1 };               % Coefficient used in right hand side.

  % Define boundary conditions.
  if( any(strcmp(opt.sfun,{'sf_tri_H3','sf_quad_H3'})) )   % Prescribed derivatives at end points for Hermite elements.
    fea.bdr.d = {{ 0 0 0 0 ;
                   0 '-ex_poisson3_derivative(y)' 0 'ex_poisson3_derivative(y)' ;
                   'ex_poisson3_derivative(x)' 0 '-ex_poisson3_derivative(x)' 0 }};
  else
    fea.bdr.d     = cell(1,n_bdr);
   [fea.bdr.d{:}] = deal(0);          % Assign zero to all boundaries (Dirichlet).
  end
  fea.bdr.n     = cell(1,n_bdr);    % No Neumann boundaries ('fea.bdr.n' empty).

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid,'icub',opt.icub);   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  g = rectgrid( 20 );
  u = evalexpr( 'u', g.p, fea );
  fv.faces    = g.c';
  fv.vertices = [g.p' u];
  fv.facevertexcdata = u;
  fv.facecolor = 'interp';
  patch( fv )
  grid on, axis normal, view(3)
  xlabel( 'x' )
  ylabel( 'y' )
  title('Solution u')
end


% Error checking.
x = linspace( 0, 1, 11 );
[x,y] = meshgrid(x,x);
u = evalexpr( 'u', [x(:) y(:)]', fea );
u_ref = l_poisol2( x(:), y(:), 6 );
u_diff = u - u_ref;
err = norm(u_diff);
if( ~isempty(fid) )
  fprintf(fid,'\nL2  Error: %f\n',err)
  fprintf(fid,'\nL00 Error: %f\n',max(abs(u_diff)))
  fprintf(fid,'\n\n')
end


out.err  = err;
out.tol  = opt.tol;
out.pass = out.err<opt.tol;
if ( nargout==0 )
  clear fea out
end


%   -----------------------------------------------------------------------------
function [u,n]=l_poisol2(x,y,n_dacc)
% Solution to Poisson equation on a unit square with n_dacc digits accuracy.

xr = linspace( 0, 1, 11 );
[xr,yr] = meshgrid(xr,xr);
if( isequal(x,xr(:)) && isequal(y,yr(:)) )
  u = refsol();
  n = inf;
  return
end

n_max = 1000;

u = zeros(size(x));
for n=1:2:n_max

  u0 = u;

  n2    = n^2;
  const = (2/pi)^4/n;
  sinj  = sin(n*pi*y);
  sini  = sin(n*pi*x);
  for i=1:2:n
    fac = const/(i*(i^2+n2));
    u   = u + fac*( sinj.*sin(i*pi*x) + sini.*sin(i*pi*y) );
  end
  u = u - const/(2*n^3)*sini.*sinj;

  udiff = floor(u0*10^n_dacc)/(10^n_dacc) - ...
          floor( u*10^n_dacc)/(10^n_dacc);
  if ( any(udiff(:)) )
    continue
  end

  break

end

%   -----------------------------------------------------------------------------
function [ u ] = refsol()
u = [ 0;
      0;
      0;
      0;
      0;
      0;
      0;
      0;
      0;
      0;
      0;
      0;
      1.3071453e-002;
      2.0880607e-002;
      2.5627979e-002;
      2.8217102e-002;
      2.9042039e-002;
      2.8217102e-002;
      2.5627979e-002;
      2.0880607e-002;
      1.3071453e-002;
      6.8917302e-018;
      0;
      2.0880607e-002;
      3.4646992e-002;
      4.3341187e-002;
      4.8155572e-002;
      4.9698164e-002;
      4.8155572e-002;
      4.3341187e-002;
      3.4646992e-002;
      2.0880607e-002;
      1.0121302e-017;
      0;
      2.5627979e-002;
      4.3341187e-002;
      5.4841059e-002;
      6.1298687e-002;
      6.3379719e-002;
      6.1298687e-002;
      5.4841059e-002;
      4.3341187e-002;
      2.5627979e-002;
      1.2034096e-017;
      0;
      2.8217102e-002;
      4.8155572e-002;
      6.1298687e-002;
      6.8744311e-002;
      7.1153117e-002;
      6.8744311e-002;
      6.1298687e-002;
      4.8155572e-002;
      2.8217102e-002;
      1.3069633e-017;
      0;
      2.9042039e-002;
      4.9698164e-002;
      6.3379719e-002;
      7.1153117e-002;
      7.3671353e-002;
      7.1153117e-002;
      6.3379719e-002;
      4.9698164e-002;
      2.9042039e-002;
      1.3398766e-017;
      0;
      2.8217102e-002;
      4.8155572e-002;
      6.1298687e-002;
      6.8744311e-002;
      7.1153117e-002;
      6.8744311e-002;
      6.1298687e-002;
      4.8155572e-002;
      2.8217102e-002;
      1.3069633e-017;
      0;
      2.5627979e-002;
      4.3341187e-002;
      5.4841059e-002;
      6.1298687e-002;
      6.3379719e-002;
      6.1298687e-002;
      5.4841059e-002;
      4.3341187e-002;
      2.5627979e-002;
      1.2034096e-017;
      0;
      2.0880607e-002;
      3.4646992e-002;
      4.3341187e-002;
      4.8155572e-002;
      4.9698164e-002;
      4.8155572e-002;
      4.3341187e-002;
      3.4646992e-002;
      2.0880607e-002;
      1.0121302e-017;
      0;
      1.3071453e-002;
      2.0880607e-002;
      2.5627979e-002;
      2.8217102e-002;
      2.9042039e-002;
      2.8217102e-002;
      2.5627979e-002;
      2.0880607e-002;
      1.3071453e-002;
      6.8917302e-018;
      0;
      6.8917302e-018;
      1.0121302e-017;
      1.2034096e-017;
      1.3069633e-017;
      1.3398766e-017;
      1.3069633e-017;
      1.2034096e-017;
      1.0121302e-017;
      6.8917302e-018;
      1.0896965e-032];
