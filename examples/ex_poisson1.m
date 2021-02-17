function [ fea, out ] = ex_poisson1( varargin )
%EX_POISSON1  1D Poisson equation example.
%
%   [ FEA, OUT ] = EX_POISSON1( VARARGIN ) Poisson equation on a line with a
%   constant source term equal to 1, homogenous boundary conditions, and
%   exact solution (-x^2+x)/2. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {1/10}          Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       iphys       scalar 0/{1}           Use physics mode to define problem    (=1)
%                                          or directly define fea.eqn/bdr fields (=0)
%                                          or use core assembly functions        (<0)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',   1/10;
            'sfun',   'sflag1';
            'refsol', '(-x^2+x)/2';
            'fsrc',   '1';
            'iphys',  1;
            'icub',   2;
            'iplot',  1;
            'tol',    2e-2;
            'fid',    1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Grid generation.
if( opt.hmax>0 )
  nx = round( 1/opt.hmax );
  fea.grid = linegrid( nx, 0, 1 );
else   % Scrambled testing grid.
  fea.grid.p = [ 0 1/10 4/10 1/3 1 1-1/3 ];
  fea.grid.c = [ 1 4 2 6 3 ;
                 2 3 4 5 6 ];
  fea.grid.a = [ 0 2 1 4 3 ;
                 2 4 3 0 5 ];
  fea.grid.b = [ 1 1 1 -1 ;
                 4 2 2  1 ]';
  fea.grid.s = ones(1,5);
end
n_bdr = 2;   % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' };                    % Coordinate name.
switch opt.iphys

 case 0   % Directly define fea.eqn/bdr fields.

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.

  % Define equation system.
  fea.eqn.a.form = { [2;2] };           % First row indicates test function space   (2=x-derivative),
                                        % second row indicates trial function space (2=x-derivative).
  fea.eqn.a.coef = { 1 };               % Coefficient used in assembling stiffness matrix.

  fea.eqn.f.form = { 1 };               % Test function space to evaluate in right hand side (1=function values).
  fea.eqn.f.coef = { opt.fsrc };        % Coefficient used in right hand side.

  % Define boundary conditions.
  if( strcmp(opt.sfun(end-1:end),'H3') )   % Prescribed derivatives at end points for Hermite elements.
    fea.bdr.d = {{ 0 0 ; 1/2 -1/2 }};
  else
    fea.bdr.d     = cell(1,n_bdr);
    [fea.bdr.d{:}] = deal(0);              % Assign zero to all boundaries (homogenous Dirichlet conditions).
  end

  fea.bdr.n     = cell(1,n_bdr);        % No Neumann boundaries ('fea.bdr.n' empty).

  % Parse and solve problem.
  fea       = parseprob(fea);           % Check and parse problem struct.
  fea.sol.u = solvestat(fea,'fid',fid,'icub',opt.icub); % Call to stationary solver.

 case 1   % Use physics mode.

  fea = addphys(fea,@poisson);          % Add Poisson equation physics mode.
  fea.phys.poi.sfun = { opt.sfun };     % Set shape function.
  fea.phys.poi.eqn.coef{3,4} = { opt.fsrc };   % Set source term coefficient.
  fea.phys.poi.bdr.coef{1,end} = repmat({0},1,n_bdr);   % Set Dirichlet boundary coefficient to zero.
  fea = parsephys(fea);                 % Check and parse physics modes.
  if( strcmp(opt.sfun(end-1:end),'H3') )   % Prescribed derivatives at end points for Hermite elements.
    fea.bdr.d = {{ 0 0 ; 1/2 -1/2 }};
  end

  % Parse and solve problem.
  fea       = parseprob(fea);           % Check and parse problem struct.
  fea.sol.u = solvestat(fea,'fid',fid,'icub',opt.icub); % Call to stationary solver.

 otherwise   % Use core assembly functions.

  fea.dvar  = { 'u' };                  % Dependent variable name.
  fea.sfun  = { opt.sfun  };            % Shape function.
  fea       = parseprob(fea);           % Check and parse problem struct.

  % Assemble stiffness matrix.
  form  = [2;2];
  sfun  = {opt.sfun;opt.sfun};
  coefa = 1;
  sind  = 1;
  i_cub = opt.icub;

  [vRowInds,vColInds,vAvals,n_rows,n_cols] = ...
    assemblea(form,sfun,coefa,i_cub,fea.grid.p,fea.grid.c,fea.grid.a,fea.grid.s,[]);
  A = sparse(vRowInds,vColInds,vAvals,n_rows,n_cols);

  % Check and compare with finite difference stencil.
  if (strcmp(opt.sfun,'sflag1'))
    h = 1/nx;
    n = nx+1;
    e = ones(n,1);
    A_ref = 1/h*spdiags([-e 2*e -e], -1:1, n, n);
    A_ref(1) = A_ref(1)/2;
    A_ref(end) = A_ref(end)/2;

    err = norm(A(:)-A_ref(:));
    if err>opt.tol
      out.err  = err;
      out.pass = -1;
      return
    end
  end

  form  = 1;
  sfun  = sfun{1};
  coeff = 1;

  f = assemblef(form,sfun,coeff,i_cub,fea.grid.p,fea.grid.c,fea.grid.a,fea.grid.s,[]);


  % Check and compare with finite difference stencil.
  if (strcmp(opt.sfun,'sflag1'))
    f_ref          = coeff*h*ones(n,1);
    f_ref([1 end]) = coeff*1/2*h;

    err = norm(f-f_ref);
    if err>1e-6
      out.err  = err;
      out.pass = -2;
      return
    end
  end

  % Set homogenous Dirichlet boundary conditions on first and last dof/node.
  bind = [1 nx+1];
  A = A';   %'
  A(:,bind) = 0;         % Zero out Dirichlet BC rows.
  for i=1:length(bind)   % Loop to set diagonal entry to 1.
    i_a = bind(i);
    A(i_a,i_a) = 1;
  end
  A = A';   %'
  f(bind) = 0;           % Set corresponding source term entries to Dirichlet BC values.

  % Solve problem.
  fea.sol.u = A\f;

end


% Postprocessing.
if ( opt.iplot>0 )
  x = linspace( 0, 1, 41 );
  u = evalexpr( 'u', x, fea )';
  figure
  subplot(3,1,1)
  plot( x, u )
  axis( [0 1 0 0.2])
  grid on
  title('Solution u')
  subplot(3,1,2)
  ux = (-x.^2+x)/2;
  plot( x, ux )
  axis( [0 1 0 0.2])
  grid on
  title('Exact solution')
  subplot(3,1,3)
  plot( x, abs(ux-u) )
  title('Error')
end


% Error checking.
xi = [1/2; 1/2];
s_err = ['abs(',opt.refsol,'-u)'];
err = evalexpr0(s_err,xi,1,1:size(fea.grid.c,2),[],fea);
ref = evalexpr0('u',xi,1,1:size(fea.grid.c,2),[],fea);
err = sqrt(sum(err.^2)/sum(ref.^2));

if( ~isempty(fid) )
  fprintf(fid,'\nL2 Error: %e\n',err)
  fprintf(fid,'\n\n')
end

out.err  = err;
out.tol  = opt.tol;
out.pass = out.err<out.tol;
if ( nargout==0 )
  clear fea out
end
