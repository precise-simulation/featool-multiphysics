function [ fea, out ] = ex_multiphase2( varargin )
%EX_MULTIPHASE2 Multiphase flow rising bubble example.
%
%   [ FEA, OUT ] = EX_MULTIPHASE2( VARARGIN ) A benchmark multiphase
%   et. al Quantitative benchmark computations of two-dimensional
%   bubble dynamics. Int. J. Numer. Meth. Fluids, 60: 1259–1288,
%   2009. doi: 10.1002/fld.1934. Accepts the following property/value
%   pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho1        scalar {1e2}           Density of the bubble
%       rho2        scalar {1e3}           Density of the surrounding fluid
%       miu1        scalar {1}             Viscosity of the bubble
%       miu2        scalar {1e1}           Viscosity of the surrounding fluid
%       gy          scalar {0.98}          Gravitational constant
%       sigma       scalar {24.5}          Coefficient of surface tension
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {1/20}          Max grid cell size
%       dt          scalar {0.1}           Time step size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       sf_c        string {sflag1}        Shape function for the level set function
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'rho1',     1e2; ...
  'rho2',     1e3; ...
  'miu1',     1; ...
  'miu2',     1e1; ...
  'gy',       0.98; ...
  'sigma',    24.5; ...
  'igrid',    0; ...
  'hmax',     1/20; ...
  'dt',       0.1; ...
  'sf_u',     'sflag2'; ...
  'sf_p',     'sflag1'; ...
  'sf_c',     'sflag1'; ...
  'iplot',    1; ...
  'itest',    0; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };   % Coordinate names.
fea.grid = rectgrid(round(1/opt.hmax),round(2/opt.hmax),[0 1;0 2]);
if ( opt.igrid==1 )
  fea.grid = quad2tri( fea.grid );
end
n_bdr = max(fea.grid.b(3,:)) + 1;   % Increment number of boundaries.
fea.grid.b(3,1) = n_bdr;            % Set first boundary edge to boundary n_bdr.


% Problem definition.
dh    = num2str(1.5*opt.hmax);      % Smoothing region width.
dw    = ['(c/',dh,'/(sqrt(cx^2+cy^2+eps)))'];
nx    = ['(cx/(sqrt(cx^2+cy^2+eps)))'];
ny    = ['(cy/(sqrt(cx^2+cy^2+eps)))'];
smhs  = ['((0.5*(1+',dw,'+1/pi*sin(pi*',dw,')))*(',dw,'>-1)*(',dw,'<1)+(',dw,'>=1))'];   % Smooth heaviside function.
smdel = ['0.5*(1+cos(pi*',dw,'))/',dh,'*(',dw,'>-1)*(',dw,'<1)'];                        % Smooth delta function.
rho   = [num2str(opt.rho1),'+',num2str(opt.rho2-opt.rho1),'*',smhs];
miu   = [num2str(opt.miu1),'+',num2str(opt.miu2-opt.miu1),'*',smhs];
sigma = num2str(opt.sigma);

% Add Navier-Stokes equations physics mode.
fea   = addphys(fea,@navierstokes);
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
fea.phys.ns.bdr.sel(n_bdr)  = 4;   % Set pressure to zero on last boundary segment.
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };

% Source terms for gravity and surface tension effects.
fg    = ['-(',rho,')*',num2str(opt.gy)];
fx1   = ['-',sigma,'*(1-(',nx,')^2)*',smdel];
fx2   = ['-',sigma,'*(-(',nx,')*(',ny,'))*',smdel];
fy1   = ['-',sigma,'*(-(',ny,')*(',nx,'))*',smdel];
fy2   = ['-',sigma,'*(1-(',ny,')^2)*',smdel];

% Add convection and diffusion physics mode for the level set equation.
fea = addphys(fea,@convectiondiffusion);
fea.phys.cd.sfun          = { opt.sf_c };
fea.phys.cd.eqn.coef{2,4} = { 0.001 };    % Add stabilizing diffusion.
fea.phys.cd.eqn.coef{3,4} = { 'u' };      % Convection velocity in the x-direction.
fea.phys.cd.eqn.coef{4,4} = { 'v' };      % Convection velocity in the y-direction.


% Parse physics modes.
fea = parsephys(fea);

% Correct source terms.
fea.eqn.f.form{1} = [     2   3];
fea.eqn.f.form{2} = [ 1   2   3];
fea.eqn.f.coef{1} = {   fx1 fx2};
fea.eqn.f.coef{2} = {fg fy1 fy2};

% Implement slip boundary conditions on vertical walls.
fea.bdr.d{2}{2} = [];
fea.bdr.d{2}{4} = [];
fea.bdr.n{2}{2} = 0;
fea.bdr.n{2}{4} = 0;

% Parse problem.
fea = parseprob(fea);


% Call to time-dependent solver.
init = { '0', '0', '0', 'sqrt((x-0.5)^2+(y-0.5)^2)-0.25' };
fea.sol.u = solvetime( fea, ...
                       'fid',     fid, ...
                       'tmax',    3.0*(~opt.itest), ...
                       'tstep',   opt.dt, ...
                       'maxnit',  24*(~opt.itest)+1, ...
                       'icub',    3, ...
                       'init',    init, ...
                       'ischeme', 3 );

% Postprocessing.
if ( opt.iplot>0 )
  subplot(2,2,1)
  postplot(fea,'surfexpr','sqrt(u^2+v^2)')
  title('Velocity field')
  subplot(2,2,2)
  postplot(fea,'surfexpr','p','evalstyle','exaxct')
  title('Pressure')
  subplot(2,2,3)
  postplot(fea,'surfexpr','c')
  title('Level set field')
  subplot(2,2,4)
  postplot(fea,'isoexpr','c','isolev',[0 0],'arrowexpr',{'u','v'})
  plot([0 1 1 0 0],[0 0 2 2 0],'k')
  title('Interface')
end


out.err  = [];
out.pass = [];
if ( nargout==0 )
  clear fea out
end
