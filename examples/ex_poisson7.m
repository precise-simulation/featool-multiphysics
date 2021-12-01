function [ fea, out ] = ex_poisson7( varargin )
%EX_POISSON7  Poisson equation on a unit circle with a point source.
%
%   [ FEA, OUT ] = EX_POISSON7( VARARGIN ) Poisson equation on a unit circle with
%   a point source (represented by a point constraint) and exact solution
%   u = -1/(2*pi)*log(r). Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.1}           Max grid cell size (<0 quadrilateral grid)
%       sfun        string {sflag1}        Shape function
%       iphys       scalar 0/{1}           Use physics mode to define problem    (=1)
%                                          or directly define fea.eqn/bdr fields (=0)
%       solver      string fenics/{default} Use FEniCS or default solver
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'hmax',     0.02;
  'refsol',   '-1/(2*pi)*log(sqrt(x^2+y^2))';
  'sfun',     'sflag1';
  'iphys',    1;
  'solver',   '';
  'iplot',    1;
  'tol',      0.2;
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.geom.objects = { gobj_circle() };


% Grid generation.
if( opt.hmax<0 )
  fea.grid = circgrid(10);
else
  hmax      = opt.hmax;
  fh        = @(p,varargin) 3*hmax + (p(:,1).^2+p(:,2).^2);
  fea.grid  = gridgen( fea, 'hmax', hmax, 'fid', fid, 'fixpnt', [0 0], 'hdfcn', fh ) ;
end


% Problem definition.
fea.sdim  = { 'x' 'y' };                % Coordinate names.
if( opt.iphys==1 )

  fea = addphys( fea, @poisson );       % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { 0 };   % Set source term coefficient.
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
  n_bdr         = max(fea.grid.b(3,:));
  fea.bdr.d     = cell(1,n_bdr);
 [fea.bdr.d{:}] = deal(0);              % Assign zero to all boundaries (Dirichlet).

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

end


% Set point constraint.
[~,i_mid] = min( fea.grid.p(1,:).^2 + fea.grid.p(2,:).^2 );
fea.pnt.index = i_mid;
fea.pnt.type  = 'source';
fea.pnt.dvar  = 1;
fea.pnt.expr  = 1;


% Parse and solve problem.
fea = parseprob( fea );               % Check and parse problem struct.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea );
else
  fea.sol.u = solvestat( fea, 'fid', fid );   % Call to stationary solver.
end


% Postprocessing.
s_err = ['abs(',opt.refsol,'-u)'];
if( opt.iplot>0 )
  figure
  subplot(3,1,1)
  postplot( fea, 'surfexpr', 'u', 'surfhexpr', 'u', 'axequal', 'on' )
  title('Solution u')

  subplot(3,1,2)
  postplot( fea, 'surfexpr', opt.refsol, 'surfhexpr', opt.refsol, 'axequal', 'on' )
  title('Exact solution')

  subplot(3,1,3)
  postplot( fea, 'surfexpr', s_err, 'surfhexpr', s_err, 'axequal', 'on' )
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
out.pass = out.err<opt.tol;
if( nargout==0 )
  clear fea out
end
