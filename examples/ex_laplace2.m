function [ fea, out ] = ex_laplace2( varargin )
%EX_LAPLACE2 2D Laplace equation on a circle with nonzero boundary conditions.
%
%   [ FEA, OUT ] = EX_LAPLACE2( VARARGIN ) Laplace equation on a circle with
%   source term 0, nonzero boundary conditions and exact solution u=r^3*sin(3*th).
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       radi        scalar {1}             Radius of circle
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
  'radi',     1; ...
  'hmax',     0.1; ...
  'refsol',   'hypot(x,y)^3*sin(3*atan2(y,x))'; ...
  'sfun',     'sflag1'; ...
  'icub',     2; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_circle( [0 0], opt.radi );
fea.geom.objects = { gobj };


% Grid generation.
fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid,'dprim',false);
n_bdr    = max(fea.grid.b(3,:));        % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };                % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { 0 };   % Set source term to zero.
  fea.phys.poi.bdr.coef{1,end} = repmat({opt.refsol},1,n_bdr);
  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2 3;2 3] };       % First row indicates test function space   (2=x-derivative + 3=y-derivative),
                                        % second row indicates trial function space (2=x-derivative + 3=y-derivative).
  fea.eqn.a.coef = { 1 };               % Coefficient used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { 0 };               % Coefficient used in right hand side.

  % Define boundary conditions.
  fea.bdr.d     = cell(1,n_bdr);
 [fea.bdr.d{:}] = deal(opt.refsol);     % Assign reference solution to all boundaries (Dirichlet).

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid,'icub',opt.icub);   % Call to stationary solver.


% Postprocessing.
s_err = ['abs(',opt.refsol,'-u)'];
if ( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot(fea,'surfexpr','u','axequal','on')
  title('Solution u')
  subplot(3,1,2)
  postplot(fea,'surfexpr',opt.refsol,'axequal','on')
  title('Exact solution')
  subplot(3,1,3)
  postplot(fea,'surfexpr',s_err,'axequal','on')
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
out.pass = out.err<5e-3;
if ( nargout==0 )
  clear fea out
end
