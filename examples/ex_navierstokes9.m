function [ fea, out ] = ex_navierstokes9( varargin )
%EX_NAVIERSTOKES9 2D Axisymmetric jet impingement with heat transfer.
%
%   [ FEA, OUT ] = EX_NAVIERSTOKES9( VARARGIN ) Sets up and solves stationary axisymmetric
%   flow of a jet impacting in a thin constrained region with a coupled heat transfer mode.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       Re          scalar {100}           Jet Reynolds number
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.2}           Max grid cell size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'Re',       100; ...
  'igrid',    0; ...
  'hmax',     0.2; ...
  'sf_u',     'sflag1'; ...
  'sf_p',     'sflag1'; ...
  'iplot',    1; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
r = 10;   % Radius.
l = 3;    % Length.
gobj = gobj_polygon( [0 0; r 0; r l; 1 l; 0 l] );
fea.geom.objects = { gobj };
fea.sdim = { 'r' 'z' };   % Space coordinate names.


% Grid generation.
if ( opt.igrid==1 )
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
else
  fea.grid = rectgrid( round(r/opt.hmax), round(l/opt.hmax), [0 r;0 l] );
  if( opt.igrid<0 )
    fea.grid = quad2tri( fea.grid );
  end
  ic = fea.grid.b(1,:);
  ii = fea.grid.b(2,:);
  jj = mod(ii,size(fea.grid.c,1)) + 1;
  p1 = fea.grid.p( :, fea.grid.c(sub2ind(size(fea.grid.c),ii,ic)) );
  p2 = fea.grid.p( :, fea.grid.c(sub2ind(size(fea.grid.c),jj,ic)) );
  pm = [ p1 + p2 ]'/2;

  ix4 = find( (pm(:,1)<1) .* (abs(pm(:,2)-3)<eps) );
  fea.grid.b(3,ix4) = 4;

  ix5 = find( pm(:,1)<=eps );
  fea.grid.b(3,ix5) = 5;
end


% Boundary specifications.
i_plate    = 1;
i_inflow   = 4;
i_outflow  = 2;
i_symmetry = 5;
inflow_bc  = -opt.Re;


% Problem definition.

% Add Navier-Stokes equations physics mode.
fea = addphys( fea, {@navierstokes 1} );
fea.phys.ns.eqn.coef{1,end} = { 1 };
fea.phys.ns.eqn.coef{2,end} = { 1 };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };


fea.phys.ns.bdr.sel(i_inflow)   = 2;
fea.phys.ns.bdr.sel(i_outflow)  = 3;
fea.phys.ns.bdr.coef{2,end}{2,i_inflow} = inflow_bc;
fea.phys.ns.bdr.sel(i_symmetry) = 5;

% Add heat transfer physics mode.
fea = addphys( fea, {@heattransfer 1} );
fea.phys.ht.eqn.coef{4,end} = fea.phys.ns.dvar{1};
fea.phys.ht.eqn.coef{5,end} = fea.phys.ns.dvar{2};
fea.phys.ht.sfun            = { opt.sf_u };

fea.phys.ht.bdr.sel(:)         = 3;
fea.phys.ht.bdr.sel(i_inflow)  = 1;
fea.phys.ht.bdr.sel(i_plate)   = 1;
fea.phys.ht.bdr.sel(i_outflow) = 2;
fea.phys.ht.bdr.coef{1,end}{1,i_inflow} = 1;

% Parse physics modes.
fea = parsephys(fea);


% Parse and solve problem.
fea       = parseprob( fea );   % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'maxnit', 50, 'fid', fid );


% Error checking.
u = evalexpr( 'u', fea.grid.p, fea );
w = evalexpr( 'w', fea.grid.p, fea );
out.err(1) = (max(u) - 75.9)/75.9;
out.err(2) = (min(u) +  9.8)/(-9.8);
out.err(3) = (max(w) -  7.9)/7.9;
out.err(4) = (min(w) +  106)/(-106);
out.pass   = all( out.err<0.1 );


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+w^2)', ...
                 'arrowexpr', {'u' 'w'}, 'arrowspacing', [60 20], ...
                 'isoexpr', 'sqrt(u^2+w^2)' )
  title( 'Velocity field' )

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'T', 'isoexpr', 'T' )
  title( 'Temperature' )
end


if( nargout==0 )
  clear fea out
end
