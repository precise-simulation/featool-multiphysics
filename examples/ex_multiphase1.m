function [ fea, out ] = ex_multiphase1( varargin )
%EX_MULTIPHASE1 Static bubble example.
%
%   [ FEA, OUT ] = EX_MULTIPHASE1( VARARGIN ) Sets up and solves a stationary
%   static bubble example where the pressure jump should be equal to
%   sigma/r. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho         scalar {1e4}           Density
%       miu         scalar {1}             Molecular/dynamic viscosity
%       sigma       scalar {1}             Coefficient of surface tension
%       r           scalar {0.25}          Bubble radius
%       lbi         scalar 1/{0}           Exact or transformed surface tension source term
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {1/30}          Max grid cell size
%       ischeme     scalar {3}             Time stepping scheme
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
  'rho',      1e4; ...
  'miu',      1; ...
  'sigma',    1; ...
  'r',        0.25; ...
  'lbi',      0; ...
  'igrid',    1; ...
  'hmax',     1/30; ...
  'ischeme'   3; ...
  'sf_u',     'sflag2'; ...
  'sf_p',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      0.3; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };   % Coordinate names.
fea.grid = rectgrid(round(1/opt.hmax),round(2/opt.hmax),[0 1;0 1]);
if ( opt.igrid==1 )
  fea.grid = quad2tri( fea.grid );
end
n_bdr = max(fea.grid.b(3,:)) + 1;   % Increment number of boundaries.
fea.grid.b(3,1) = n_bdr;            % Set first boundary edge to boundary n_bdr.


% Problem definition.
phi   = ['(sqrt((x-0.5)^2+(y-0.5)^2)-',num2str(opt.r),')'];   % Level set function.
dh    = num2str(1*opt.hmax);   % Smoothing region width.
dw    = ['(',phi,'/',dh,')'];
smhs  = ['((0.5*(1+',dw,'+1/pi*sin(pi*',dw,')))*(',dw,'>-1)*(',dw,'<1)+(',dw,'>=1))'];   % Smooth heaviside function.
smdel = ['0.5*(1+cos(pi*',dw,'))/',dh,'*(',dw,'>-1)*(',dw,'<1)'];                        % Smooth delta function.
sigma = num2str(opt.sigma);
nx    = ['(2*x-1)/(2*sqrt((x-0.5)^2+(y-0.5)^2+eps))'];
ny    = ['(2*y-1)/(2*sqrt((x-0.5)^2+(y-0.5)^2+eps))'];

% Surface tension force source term.
if ( ~opt.lbi )   % Exact curvature.
  kappa = ['1/sqrt((x-0.5)^2+(y-0.5)^2+eps)'];
  fx1   = ['-',sigma,'*',kappa,'*',nx,'*',smdel];
  fy1   = ['-',sigma,'*',kappa,'*',ny,'*',smdel];
else   % With Laplace-Beltrami transformation.
  fx1   = ['-',sigma,'*(1-(',nx,')^2)*',smdel];
  fx2   = ['-',sigma,'*(-(',nx,')*(',ny,'))*',smdel];
  fy1   = ['-',sigma,'*(-(',ny,')*(',nx,'))*',smdel];
  fy2   = ['-',sigma,'*(1-(',ny,')^2)*',smdel];
end


% Add Navier-Stokes equations physics mode.
fea = addphys(fea,@navierstokes);
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };
fea.phys.ns.eqn.coef{3,end} = { fx1 };
fea.phys.ns.eqn.coef{4,end} = { fy1 };
fea.phys.ns.bdr.sel(n_bdr)  = 4;   % Set pressure to zero on last boundary segment.


% Parse problem.
fea = parsephys(fea);
if ( opt.lbi )
  fea.eqn.f.form{1} = [2 3];
  fea.eqn.f.form{2} = [2 3];
  fea.eqn.f.coef{1} = {fx1 fx2};
  fea.eqn.f.coef{2} = {fy1 fy2};
end
fea = parseprob(fea);


% Call to time-dependent solver.
fea.sol.u = solvetime( fea, ...
                       'fid',     fid, ...
                       'tmax',    0.3, ...
                       'tstep',   0.1, ...
                       'icub',    3, ...
                       'ischeme', opt.ischeme );

% Postprocessing.
if ( opt.iplot>0 )
  postplot(fea,'surfexpr','p')
  title('Pressure')
end


% Error checking.
ind_p    = [fea.eqn.ndof(1)+fea.eqn.ndof(2)+1:fea.eqn.ndof(1)+fea.eqn.ndof(2)+fea.eqn.ndof(3)];
dp       = max(fea.sol.u(ind_p,end)) - min(fea.sol.u(ind_p,end));
out.err  = abs(opt.sigma/opt.r - dp)/(opt.sigma/opt.r);
out.pass = out.err<opt.tol;

if ( nargout==0 )
  clear fea out
end
