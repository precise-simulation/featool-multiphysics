function [ fea, out ] = ex_planestress7( varargin )
%EX_PLANESTRESS7 Thermally induced displacements on a rectangle
%
%   [ FEA, OUT ] = EX_PLANESTRESS7( VARARGIN ) Thermally induced
%   displacements on a rectangle
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       pss         logical {true}         Use plane stress (or plane strain)
%       hmax        scalar {1/20}          Max grid cell size
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, ~0=triangles)
%       sfun        string {sflag1}        Shape function for displacements
%       solver      string {}              Solver selection default, fenics
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'pss',    true;
            'hmax',   1/20;
            'igrid',  1;
            'sfun',   'sflag1';
            'solver', '';
            'iplot',  1;
            'tol',    1e-4;
            'fid',    1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry and grid.
fea.sdim = { 'x' 'y' };   % Coordinate names.
fea.grid = rectgrid(1/opt.hmax, 1/opt.hmax, [0 0.01; 0 0.01]);
if( opt.igrid~=0 || strcmpi(opt.solver,'fenics') )
  fea.grid = quad2tri( fea.grid );
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Problem definition.
if( opt.pss )
  fea = addphys( fea, @planestress );
  mode = 'pss';
else
  fea = addphys( fea, @planestrain );
  mode = 'psn';
end
fea.phys.(mode).eqn.coef{1,end} = { 0.3 };
fea.phys.(mode).eqn.coef{2,end} = { 100e9 };
fea.phys.(mode).eqn.coef{6,end} = { 1e-3 };
fea.phys.(mode).eqn.coef{7,end} = { 1 };
fea.phys.(mode).sfun            = { opt.sfun opt.sfun };


% Boundary conditions.
dtol = 1e-6;
lbdr = findbdr( fea, ['x<',num2str(dtol)] );
bbdr = findbdr( fea, ['y<',num2str(dtol)] );
n_bdr  = max(fea.grid.b(3,:));
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
[bctype{1,lbdr}] = deal(1);
[bctype{2,bbdr}] = deal(1);
fea.phys.(mode).bdr.coef{1,end} = bccoef;
fea.phys.(mode).bdr.coef{1,5}   = bctype;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );                          % Check and parse problem struct.
if( strcmp(opt.solver,'fenics') )
  fea = fenics(fea,'fid',fid);
else
  fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.
end


% Postprocessing.
s_disp = fea.phys.(mode).eqn.vars{2,end};
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_disp )
  title( 'Total displacement' )
end


% Error checking.
if( opt.pss )
  u_max_ref = 0;
  u_mag_max_ref = 1.41421e-5;
else
  u_max_ref = 1.3e-5;
  u_mag_max_ref = 1.83848e-5;
end
[u_min,u_max] = minmaxsubd('u',fea);
[v_min,v_max] = minmaxsubd('v',fea);
[u_mag_min,u_mag_max] = minmaxsubd('sqrt(u^2+v^2)',fea);
out.err = [ abs(u_min), ...
            abs(u_max-u_max_ref), ...
            abs(v_min), ...
            abs(v_max-u_max_ref), ...
            abs(u_mag_min), ...
            abs(u_mag_max-u_mag_max_ref)/u_mag_max_ref ];
out.pass = all( out.err(:) <= opt.tol );


if( nargout==0 )
  clear fea out
end
