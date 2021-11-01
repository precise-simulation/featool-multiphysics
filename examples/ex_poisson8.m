function [ fea, out ] = ex_poisson8( varargin )
%EX_POISSON8 2D Poisson equation example on a unit square with integral constraint.
%
%   [ FEA, OUT ] = EX_POISSON8( VARARGIN ) Poisson equation on a
%   [0..1]^2 unit square with all Neumann boundary conditions,
%   integral constraint, and exponential source term.
%
%   Ref. https://fenicsproject.org/olddocs/dolfin/dev/python/demos/neumann-poisson/demo_neumann-poisson.py.html
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Max grid cell size
%       sfun        string {sflag1}        Shape function
%       ischeme     scalar {0}             Time stepping scheme (0 = stationary)
%       solver      string fenics/{default} Use FEniCS or default solver
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
            'igrid',    0;
            'hmax',     0.1;
            'sfun',     'sflag1';
            'ischeme',  0;
            'solver',   '';
            'iplot',    1;
            'tol',      0.02;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.geom.objects = { gobj_rectangle() };


% Grid generation.
switch( opt.igrid )
  case -1
    fea.grid = rectgrid(round(1/opt.hmax),round(1/opt.hmax),[0 1;0 1]);
    fea.grid = quad2tri(fea.grid);
  case 0
    fea.grid = rectgrid(round(1/opt.hmax),round(1/opt.hmax),[0 1;0 1]);
  case 1
    fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
end
if( strcmp(opt.solver,'fenics') && size(fea.grid.c,1)==4 )
  fea.grid = quad2tri(fea.grid);
end

% Problem definition.
fea.sdim  = { 'x', 'y' };               % Coordinate names.

fea = addphys(fea,@poisson);            % Add Poisson equation physics mode.
fea.phys.poi.sfun = { opt.sfun };       % Set shape function.
fea.phys.poi.eqn.coef{3,end}{1} = '10*exp(-((x-0.5)^2+(y-0.5)^2)/0.02)';
fea.phys.poi.bdr.sel(:) = 2;
[fea.phys.poi.bdr.coef{2,end}{:}] = deal('-sin(5*x)');
fea = parsephys(fea);                   % Check and parse physics modes.

% Add integral constraint.
constr.type = 'intsubd';
constr.dvar = 'u';
constr.expr = 0;
fea.constr  = constr;


% Parse and solve problem.
fea = parseprob( fea );               % Check and parse problem struct.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', fid, 'ischeme', opt.ischeme, 'tstep', 0.1, 'tmax', 1 );
else
  if( opt.ischeme==0 )
    fea.sol.u = solvestat( fea, 'fid', fid );
  else
    fea.sol.u = solvetime( fea, 'fid', fid, 'ischeme', opt.ischeme, 'tstep', 0.1, 'tmax', 1 );
  end
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'u', 'surfhexpr', 'u' )
end


% Error checking.
[u_min,u_max] = minmaxsubd( 'u', fea );
err = mean( [ abs(u_max-u_min-1.037)/1.037, ...
              intsubd( 'u', fea ), ...
              abs(fea.sol.u(end)-1.3007)/1.3007 ] );

out.err  = err;
out.tol  = opt.tol;
out.pass = out.err<opt.tol;
if( nargout==0 )
  clear fea out
end
