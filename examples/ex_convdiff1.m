function [ fea, out ] = ex_convdiff1( varargin )
%EX_CONVDIFF1 2D Convection and diffusion equation example on a rectangle.
%
%   [ FEA, OUT ] = EX_CONVDIFF1( VARARGIN ) Convection and diffusion equation on
%   a rectangle with exact solution u_0+c1*eta+c2*(2*cd*xi+eta^2).
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {1/40}          Max grid cell size
%       a           scalar {1}             Convection velocity in x-direction
%       b           scalar {2}             Convection velocity in y-direction
%       c1          scalar {1}             Solution constant
%       c2          scalar {0.8}           Solution constant
%       cd          scalar {0.5}           Diffusion coefficient
%       sfun        string {sflag1}        Shape function
%       iphys       scalar 0/{1}           Use physics mode to define problem    (=1)
%                                          or directly define fea.eqn/bdr fields (=0)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2024 Precise Simulation, Ltd.


cOptDef = { ...
  'igrid',    0; ...
  'hmax',     1/40; ...
  'a',        1; ...
  'b',        2; ...
  'c1',       1; ...
  'c2',       0.8; ...
  'cd',       0.5; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;
xi        = [num2str(opt.a),'*x+',num2str(opt.b),'*y'];
eta       = [num2str(opt.b),'*x-',num2str(opt.a),'*y'];
refsol    = [num2str(opt.c1),'*(',eta,')+',num2str(opt.c2),'*(2*',num2str(opt.cd),'*(',xi,')+(',eta,')^2)'];


% Geometry definition.
gobj = gobj_rectangle();
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case -1
    fea.grid = rectgrid(round(1/opt.hmax));
    fea.grid = quad2tri(fea.grid);
  case 0
    fea.grid = rectgrid(round(1/opt.hmax));
  case 1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
end
n_bdr = max(fea.grid.b(3,:));           % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };                % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys(fea,@convectiondiffusion);   % Add convection and diffusion physics mode.
  fea.phys.cd.sfun = { opt.sfun };           % Set shape function.
  fea.phys.cd.eqn.coef{2,4}   = { opt.cd };  % Set diffusion coefficient.
  fea.phys.cd.eqn.coef{3,4}   = { opt.a  };  % Convection velocity in x-direction.
  fea.phys.cd.eqn.coef{4,4}   = { opt.b  };  % Convection velocity in y-direction.
  fea.phys.cd.bdr.sel         = [1 1 1 1];
  fea.phys.cd.bdr.coef{1,end} = repmat({refsol},1,n_bdr);   % Set Dirichlet boundary coefficient to reference solution.
  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'c' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2 3 2 3;2 3 1 1] };   % First row indicates test function space   (2=x-derivative + 3=y-derivative),
                                            % second row indicates trial function space (2=x-derivative + 3=y-derivative).
  fea.eqn.a.coef = { [opt.cd opt.cd opt.a opt.b] };   % Coefficients used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { 0 };               % Coefficient used in right hand side.

  % Define boundary conditions.
  fea.bdr.d     = cell(1,n_bdr);
 [fea.bdr.d{:}] = deal(refsol);         % Assign reference solution to all boundaries (Dirichlet).

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
s_err = ['abs(',refsol,'-c)'];
if ( opt.iplot>0 )
  figure
  subplot(2,1,1)
  postplot(fea,'surfexpr','c','isoexpr','c')
  title('Solution c')
  subplot(2,1,2)
  postplot(fea,'surfexpr',s_err)
  title('Error')
end


% Error checking.
if ( size(fea.grid.c,1)==4 )
  xi = [0;0];
else
  xi = [1/3;1/3;1/3];
end
err = evalexpr0(s_err,xi,1,1:size(fea.grid.c,2),[],fea);
ref = evalexpr0('c',xi,1,1:size(fea.grid.c,2),[],fea);
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
