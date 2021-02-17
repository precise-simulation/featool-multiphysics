function [ fea, out ] = ex_natural_convection( varargin )
%EX_NATURAL_CONVECTION 2D Example for natural convection of air in a square cavity.
%
%   [ FEA, OUT ] = EX_NATURAL_CONVECTION( VARARGIN ) Sets up and solves a natural convection
%   benchmark problem. Reference solutions are for example reported in G. Davis
%   "Natural convection of air in a square cavity a bench mark numerical solution",
%   IJNMF vol. 3, 249-254 (1983).
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       Ra          scalar {1e3}           Rayleigh number
%       Pr          scalar {0.71}          Prandtl number
%       l           scalar {1}             Side length of cavity
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.05}          Max grid cell size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       sf_T        string {sflag2}        Shape function for temperature
%       iphys       scalar 0/{1}           Use physics mode to define problem (=1)
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'Ra',       1e3; ...
  'Pr',       0.71; ...
  'l',        1; ...
  'igrid',    0; ...
  'hmax',     0.05; ...
  'sf_u',     'sflag2'; ...
  'sf_p',     'sflag1'; ...
  'sf_T',     'sflag2'; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'tol',      0.1; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;

% Reference values.
Ra_ref      = [ 1e3   1e4    1e5     1e6     ];
u_max_ref   = [ 3.649 16.178 34.73   64.63   ];
y_max_ref   = [ 0.813  0.823  0.855   0.850  ];
v_max_ref   = [ 3.697 19.617 68.59  219.36   ];
x_max_ref   = [ 0.178  0.119 0.066    0.0379 ];
Nu_mean_ref = [ 1.118  2.243 4.519    8.800  ];

% Model parameters.
Ra    = opt.Ra;      % Rayleigh number.
Pr    = opt.Pr;      % Prandtl number.
l     = opt.l;       % Length of rectangular domain.
sf_u  = opt.sf_u;    % FEM shape function type for velocity.
sf_p  = opt.sf_p;    % FEM shape function type for pressure.
sf_T  = opt.sf_T;    % FEM shape function type for temperature.


% Geometry definition.
fea.sdim         = { 'x' 'y' };
fea.geom.objects = { gobj_rectangle( 0, l, 0, l ) };


% Grid generation.
if( opt.igrid<=0 )
  fea.grid = rectgrid( round(l/opt.hmax), round(l/opt.hmax), [0 l;0 l] );
  if( opt.igrid<0 )
    fea.grid = quad2tri(fea.grid);
  end
else
  fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );
end


% Boundary conditions.
n_bdr = max(fea.grid.b(3,:));    % Increment number of boundaries.
dtol  = opt.hmax;
ib_l  = findbdr( fea, ['x<',num2str(dtol)] );     % Right boundary number.
ib_r  = findbdr( fea, ['x>',num2str(l-dtol)] );   % Left boundary number.
ib_b  = findbdr( fea, ['y<',num2str(dtol)] );     % Bottom boundary number.
ib_t  = findbdr( fea, ['y>',num2str(l-dtol)] );   % Top boundary number.

% Add pressure point constraint on point closest to origin.
[~,ix] = min( fea.grid.p(1,:).^2 + fea.grid.p(2,:).^2 );
fea.pnt.index = ix;
fea.pnt.type  = 'constr';
fea.pnt.dvar  = 'p';
fea.pnt.expr  = 0';


% Problem definition.
if ( opt.iphys==1 )

  fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
  fea.phys.ns.eqn.coef{2,end} = { Pr };
  fea.phys.ns.eqn.coef{4,end} = { [num2str(Ra*Pr),'*T'] };
  fea.phys.ns.sfun            = { sf_u sf_u sf_p };   % Set shape functions.

  fea = addphys(fea,@heattransfer);     % Add heat transfer physics mode.
  fea.phys.ht.sfun            = { sf_T };
  fea.phys.ht.eqn.coef{4,end} = { fea.phys.ns.dvar{1} };
  fea.phys.ht.eqn.coef{5,end} = { fea.phys.ns.dvar{2} };
  fea.phys.ht.bdr.sel([ib_l ib_r])  = 1;
  fea.phys.ht.bdr.coef{1,end}{ib_l} = 1;

  fea = parsephys(fea);                 % Check and parse physics modes.

else

  fea.dvar  = { 'u'  'v'  'p'  'T'  };       % Dependent variable name.
  fea.sfun  = { sf_u sf_u sf_p sf_T };       % Shape function.

  % Define equation system.
  cvelx = fea.dvar{1};   % Convection velocities.
  cvely = fea.dvar{2};
  fea.eqn.a.form = { [2 3 2 3;2 3 1 1]      [2;3]                   [1;2] []; ...
                     [3;2]                  [2 3 2 3;2 3 1 1]       [1;3] []; ...
                     [2;1]                  [3;1]                   []    []; ...
                     []                     []                      []    [2 3 2 3;2 3 1 1] };
  fea.eqn.a.coef = { {2*Pr Pr cvelx cvely}   Pr                     -1 []; ...
                     Pr                     {Pr 2*Pr cvelx cvely}   -1 []; ...
                     1                       1                      [] []; ...
                     []                     []                      [] {1 1 cvelx cvely} };
  fea.eqn.f.form = { 1 1 1 1 };
  fea.eqn.f.coef = { 0 [num2str(Ra*Pr),'*',fea.dvar{4}] 0 0 };

  % Define boundary conditions.
  fea.bdr.d = cell(4,n_bdr);
  [fea.bdr.d{1:2,:}] = deal(0);

  fea.bdr.d{4,ib_l}  = 1;
  fea.bdr.d{4,ib_r}  = 0;

  fea.bdr.n = cell(4,n_bdr);
  [fea.bdr.n{4,[ib_t ib_b n_bdr]}] = deal(0);
  [fea.bdr.n{3,setdiff(1:n_bdr,n_bdr)}] = deal(0);
  [fea.bdr.n{1:2,n_bdr}] = deal(0);

end


% Parse and solve problem.
fea       = parseprob(fea);             % Check and parse problem struct.
fea.sol.u = solvestat(fea,'fid',fid);   % Call to stationary solver.


% Postprocessing.
if ( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot(fea,'surfexpr','sqrt(u^2+v^2)')
  title('Velocity field')
  subplot(1,2,2)
  postplot(fea,'surfexpr','T')
  title('Temperature')
end


% Error checking.
out.err  = nan;
out.pass = nan;
iref = find( Ra==Ra_ref );
if ( ~isempty(iref) )
  n_evalution_points = 3*l/opt.hmax;
  x_eval = linspace( 0, l, n_evalution_points );
  x_mid  = l/2*ones( 1, n_evalution_points );

  u_eval = evalexpr( 'u', [x_mid; x_eval], fea );
  [u_max,ix] = max( u_eval );
  y_max = x_eval( ix );

  v_eval = evalexpr( 'v', [x_eval; x_mid], fea );
  [v_max,ix] = max( v_eval );
  x_max = x_eval( ix );

  Nu_mean = abs( intbdr( 'Tx', fea, 1, 2 ) + intbdr( 'Tx', fea, 2, 2 ) )/2;

  out.u_max   = u_max;
  out.y_max   = y_max;
  out.v_max   = v_max;
  out.x_max   = x_max;
  out.Nu_mean = Nu_mean;

  out.err = [ abs(u_max_ref(iref)-u_max)/u_max_ref(iref);
              abs(y_max_ref(iref)-y_max)/y_max_ref(iref);
              abs(v_max_ref(iref)-u_max)/v_max_ref(iref);
              abs(x_max_ref(iref)-x_max)/x_max_ref(iref);
              abs(Nu_mean_ref(iref)-Nu_mean)/Nu_mean_ref(iref) ];
  out.pass = all( out.err<opt.tol );
end

if ( nargout==0 )
  clear fea out
end
