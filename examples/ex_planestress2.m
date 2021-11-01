function [ fea, out ] = ex_planestress2( varargin )
%EX_PLANESTRESS2 NAFEMS benchmark challenge 1 plane stress example.
%
%   [ FEA, OUT ] = EX_PLANESTRESS2( VARARGIN ) NAFEMS benchmark example to calculate von Mieses stress
%   thin plate under plane stress assumption with loads all around.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {210e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {1/20}          Max grid cell size
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        210e9; ...
  'nu',       0.3; ...
  'igrid',    0; ...
  'hmax',     1/20; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj1 = gobj_rectangle( 0, 1, 0, 1, 'R1' );
fea.geom.objects = { gobj1 };
fea.sdim = { 'x' 'y' };


% Grid generation.
if( opt.igrid==1 )
  fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
else
  nx = 1/opt.hmax;
  fea.grid = rectgrid( nx, nx );
end

% Boundary conditions.
dtol = 0.1;
if( opt.igrid==1 )
  lbdr = findbdr( fea, ['x<',num2str(dtol)] );     % Left boundary number.
  rbdr = findbdr( fea, ['x>',num2str(1-dtol)] );   % Right boundary number.
  tbdr = findbdr( fea, ['y>',num2str(1-dtol)] );   % Top boundary number.
  bbdr = findbdr( fea, ['y<',num2str(dtol)] );     % Bottom boundary number.
else
  lbdr = 4;
  rbdr = 2;
  tbdr = 3;
  bbdr = 1;
end


% Add plane stress physics mode.
fea = addphys(fea,@planestress);
fea.phys.pss.eqn.coef{1,end} = { opt.nu };
fea.phys.pss.eqn.coef{2,end} = { opt.E  };
fea.phys.pss.sfun            = { opt.sfun opt.sfun };

% Set boundary condition types.
bctype = mat2cell( zeros(2,4), [1 1], [1 1 1 1] );
fea.phys.pss.bdr.coef{1,5}   = bctype;

% Set loads on boundary.
bccoef = mat2cell( zeros(2,4), [1 1], [1 1 1 1] );
bccoef{1,tbdr} = '-x';
bccoef{2,tbdr} = 'x';
bccoef{1,bbdr} = '-(1-x)';
bccoef{2,bbdr} = '1-x';
bccoef{1,lbdr} = '1-y';
bccoef{2,lbdr} = '-(1-y)';
bccoef{1,rbdr} = 'y';
bccoef{2,rbdr} = '-y';
fea.phys.pss.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea       = parsephys(fea);             % Check and parse physics modes.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
s_vm = fea.phys.pss.eqn.vars{1,2};
if ( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_vm, 'isoexpr', s_vm )
  title('von Mieses stress')
end

% Error checking.
dtol = sqrt(eps);
sm  = evalexpr( s_vm, [0.5;0.5], fea );
s01 = evalexpr( s_vm, [dtol;1],    fea );
s10 = evalexpr( s_vm, [1;dtol],    fea );
s00 = evalexpr( s_vm, [dtol;dtol], fea );
s11 = evalexpr( s_vm, [1;1],       fea );

sol = [ sm s01 s10 s00 s11 ];
ref = [  0   0   0   2   2 ];
out.err    = norm(sol-ref)/norm(ref);
out.pass   = out.err < 0.1;


if ( nargout==0 )
  clear fea out
end
