function [ fea, out ] = ex_poisson4( varargin )
%EX_POISSON4 2D Poisson equation example on a rectangle.
%
%   [ FEA, OUT ] = EX_POISSON4( VARARGIN ) Poisson equation on a rectangle with
%   source term (1+2000^2) * pi^2 * sin(pi*x) * sin(2000*pi*y) and exact solution
%   sin(pi*x) * sin(2000*pi*y). Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       lx          scalar {1e-3}          Length of domain in x-direction
%       ly          scalar {1e-3}          Length of domain in y-direction
%       hmax        scalar {ly/20}         Max grid cell size
%       refsol      string {sin(pi*x)*sin(2000*pi*y)}   Reference solution
%       fsrc        string {(1+2000^2)*pi^2*sin(pi*x)*sin(2000*pi*y)}   Source term
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
  'lx',       0.001; ...
  'ly',       0.001; ...
  'hmax',     0.001/20; ...
  'refsol',   'sin(pi*x)*sin(2000*pi*y)'; ...
  'fsrc',     '(1+2000^2)*pi^2*sin(pi*x)*sin(2000*pi*y)'; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_rectangle( 0, opt.lx, 0, opt.lx );
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case -1
    fea.grid = rectgrid(round(opt.lx/opt.hmax),round(opt.lx/opt.hmax),[0 opt.lx;0 opt.lx]);
    fea.grid = quad2tri(fea.grid);
  case 0
    fea.grid = rectgrid(round(opt.lx/opt.hmax),round(opt.lx/opt.hmax),[0 opt.lx;0 opt.lx]);
  case 1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
end
fea.grid.p(2,:) = fea.grid.p(2,:)*opt.ly/opt.lx;   % Scale y-coordinates.
n_bdr = max(fea.grid.b(3,:));           % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };                % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { opt.fsrc };   % Set source term coefficient.
  fea.phys.poi.bdr.coef{1,end} = repmat({opt.refsol},1,n_bdr);   % Set Dirichlet boundary coefficient to reference solution.
  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2 3;2 3] };       % First row indicates test function space   (2=x-derivative + 3=y-derivative),
                                        % second row indicates trial function space (2=x-derivative + 3=y-derivative).
  fea.eqn.a.coef = { 1 };               % Coefficient used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { opt.fsrc };        % Coefficient used in right hand side.

  % Define boundary conditions.
  fea.bdr.d     = cell(1,n_bdr);
 [fea.bdr.d{:}] = deal(opt.refsol);     % Assign reference solution to all boundaries (Dirichlet).

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
s_err = ['abs(',opt.refsol,'-u)'];
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr','u','axequal','off')
  title('Solution u')
  subplot(3,1,2)
  postplot(fea,'surfexpr',opt.refsol,'axequal','off')
  title('Exact solution')
  subplot(3,1,3)
  postplot(fea,'surfexpr',s_err,'axequal','off','evalstyle','exact')
  title('Error')
end


% Error checking.
if ( size(fea.grid.c,1)==4 )
  xi = [0;0];
else
  xi = [1/3;1/3;1/3];
end
err = evalexpr0(s_err,xi,1,1:size(fea.grid.c,2),[],fea);
ref = evalexpr0('u',xi,1,1:size(fea.grid.c,2),[],fea);
err = sqrt(sum(err.^2)/sum(ref.^2));

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %f\n',err)
  fprintf(fid,'\n\n')
end


out.err  = err;
out.pass = out.err<0.1;
if ( nargout==0 )
  clear fea out
end
