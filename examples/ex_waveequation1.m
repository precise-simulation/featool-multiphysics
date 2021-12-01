function [ fea, out ] = ex_waveequation1( varargin )
%EX_WAVEEQUATION1 2D Wave equation example on a circle.
%
%   [ FEA, OUT ] = EX_WAVEEQUATION1( VARARGIN ) Wave equation on a circle with
%   zero source term by variable splitting. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       radi        scalar {1}             Radius of circle
%       c2          scalar {1}             Wave speed (squared)
%       init        scalar {1-(x^2+y^2)}   Initial shape
%       tmax        scalar {1}             Stopping time
%       tstep       scalar {0.1}           Time step size
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.05}          Grid cell size
%       iexpl       scalar {0}             Use explicit (or implicit) right hand side
%       ischeme     scalar {3}             Time stepping scheme
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
  'c2',       1; ...
  'init',     '1-(x^2+y^2)'; ...
  'tmax',     1; ...
  'igrid',    1; ...
  'hmax',     0.1; ...
  'tstep',    0.05; ...
  'iexpl'     0; ...
  'ischeme'   3; ...
  'sfun',     'sflag1'; ...
  'iphys',    1; ...
  'iplot',    1; ...
  'tol',      0.125; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_circle( [0 0], opt.radi );
fea.geom.objects = { gobj };


% Grid generation.
if( opt.igrid==1 )
  fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid,'dprim',false);
else
  fea.grid = circgrid( 16, 12, opt.radi );
  if( opt.igrid<0 )
    fea.grid = quad2tri( fea.grid );
  end
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };   % Coordinate names.
if ( opt.iphys==1 )

  fea = addphys( fea, @poisson, {'u'} );  % Add Poisson equation physics mode for 'u'.
  fea = addphys( fea, @poisson, {'v'} );  % Add Poisson equation physics mode for 'v'.
  fea.phys.poi.sfun  = { opt.sfun };      % Set shape function for 'u'.
  fea.phys.poi2.sfun = { opt.sfun };      % Set shape function for 'v'.

  % Change equations.
  if( opt.iexpl==1 )
    fea.phys.poi.eqn.seqn  = ['u'' = v'];
  else
    fea.phys.poi.eqn.seqn  = ['u'' - v_t = 0'];
  end
  fea.phys.poi2.eqn.seqn = ['v'' - ',num2str(opt.c2),'*(ux_x  + uy_y) = 0'];

  % Set homogenous Dirichlet boundary coefficients.
  fea.phys.poi.bdr.coef{1,end}  = repmat({0},1,n_bdr);
  fea.phys.poi2.bdr.coef{1,end} = repmat({0},1,n_bdr);

  % Check and parse physics modes.
  fea = parsephys(fea);

else

  fea.dvar  = { 'u' 'v' };              % Dependent variable names.
  fea.sfun  = { opt.sfun opt.sfun };    % Shape functions.

  % Mass matrix.
  fea.eqn.m.form = { [1;1] []    ;
                     []    [1;1] };
  fea.eqn.m.coef = { 1  [] ;
                     [] 1  };

  % System/iteration matrix A_u := 0, A_vu = c2*( u_xx + u_yy ).
  if( opt.iexpl==1 )
    a12 = [];
  else
    a12 = [1;1];
  end
  fea.eqn.a.form = { []        a12 ;
                     [2 3;2 3] [] };
  fea.eqn.a.coef = { []     -1 ;
                     opt.c2 [] };

  % Source term  f_u := v, f_v := 0.
  fea.eqn.f.form = {  1 ; 1 };
  if( opt.iexpl==1 )
    fea.eqn.f.coef = { 'v'; 0 };
  else
    fea.eqn.f.coef = { 0 ; 0 };
  end

  % Define homogenous Dirichlet conditions.
  fea.bdr.d     = cell(2,n_bdr);
 [fea.bdr.d{:}] = deal( 0 );       % Zero Dirichlet conditions everyhere.
  fea.bdr.n     = cell(2,n_bdr);   % Clear Neumann conditions.

end


% Parse and solve problem.
fea = parseprob(fea);
[fea.sol.u,fea.sol.t] = solvetime( fea, 'fid', fid, ...
                                   'init',  { opt.init 0 }, ...
                                   'icub',    4, ...
                                   'imass',   4, ...
                                   'tmax',    opt.tmax, ...
                                   'tstep',   opt.tstep, ...
                                   'ischeme', opt.ischeme );

% Postprocessing.
xp = [0; 0];
if ( opt.iplot>0 )
  figure
  subplot(1,2,1)
  isol = numel(fea.sol.t);
  postplot( fea, 'surfexpr', 'u', 'axequal', 'on', 'solnum', isol )
  title(['Solution at time ',num2str(fea.sol.t(isol))])

  for isol=1:numel(fea.sol.t)
    u(isol) = evalexpr( 'u', xp, fea, isol );
  end
  subplot(1,2,2)
  plot( fea.sol.t, u )
  title(['Solution at point (',num2str(xp'),')'])
  ylabel('u')
  xlabel('time')
end

u_ref    = -0.958;
[~,isol] = min(abs(fea.sol.t-1));
out.err  = abs( u_ref - evalexpr( 'u', xp, fea, isol ) )/abs(u_ref);
out.pass = out.err < opt.tol;
if ( nargout==0 )
  clear fea out
end
