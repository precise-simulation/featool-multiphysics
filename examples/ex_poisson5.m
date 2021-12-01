function [ fea, out ] = ex_poisson5( varargin )
%EX_POISSON5 3D Poisson equation example on a unit sphere.
%
%   [ FEA, OUT ] = EX_POISSON5( VARARGIN ) Poisson equation on a unit sphere with
%   source term f=1 and homogenous (zero) Dirichlet boundary conditions on the
%   sphere surface. The exact solution to this problem is u_ref=(1-r^2)/6 where
%   r is the radius from the origin. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.35}          Max grid cell size
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.

cOptDef = { ...
  'igrid',    0; ...
  'hmax',     0.35; ...
  'refsol',   '(1-(x^2+y^2+z^2))/6'; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'icub',     2; ...
  'iplot',    1; ...
  'tol',      0.2;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_sphere();
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case -2
    fea.grid = blockgrid( 1, 1, 1, [-1 1;-1 1;-1 1] );
  case -1
    fea.grid = spheregrid(round(1/opt.hmax*4/7),round(1/opt.hmax*3/7));
    fea.grid = hex2tet(fea.grid);
  case 0
    fea.grid = spheregrid(max(2,5*round(1/opt.hmax*4/7)),max(1,5*round(1/opt.hmax*3/7)));
  case 1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid,'dprim',false);
end
n_bdr = max(fea.grid.b(3,:));           % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' 'z' };            % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { 1 };   % Set source term coefficient.
  fea.phys.poi.bdr.coef{1,end} = repmat({opt.refsol},1,n_bdr);   % Assign reference solution to all boundaries (Dirichlet).
  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2 3 4;2 3 4] };   % First row indicates test function space   (2=x-derivative + 3=y-derivative),
                                        % second row indicates trial function space (2=x-derivative + 3=y-derivative).
  fea.eqn.a.coef = { 1 };               % Coefficient used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { 1 };               % Coefficient used in right hand side.

  % Define boundary conditions.
  fea.bdr.d     = cell(1,n_bdr);
 [fea.bdr.d{:}] = deal(opt.refsol);     % Assign reference solution to all boundaries (Dirichlet).

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid,'icub',opt.icub);   % Call to stationary solver.


% Postprocessing.
if ( opt.iplot>0 )
  figure
  postplot(fea,'surfexpr','u','selexpr','(y>0)','axequal','on')
end


% Error checking.
s_err = ['abs(',opt.refsol,'-u)'];
if ( size(fea.grid.c,1)==8 )
  xi = [0;0;0];
else
  xi = [1/4;1/4;1/4;1/4];
end
err = evalexpr0(s_err,xi,1,1:size(fea.grid.c,2),[],fea);
ref = evalexpr0('u',xi,1,1:size(fea.grid.c,2),[],fea);
err = sqrt(sum(err.^2)/sum(ref.^2));

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.tol  = opt.tol;
out.pass = out.err<out.tol;
if ( nargout==0 )
  clear fea out
end
