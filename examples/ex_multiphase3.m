function [ fea, out ] = ex_multiphase3( varargin )
%EX_MULTIPHASE3 Breaking dam multiphase flow example.
%
%   [ FEA, OUT ] = EX_MULTIPHASE3( VARARGIN ) Breaking dam multiphase flow benchmark
%   problem as define by Martin J, Moyce W. Part IV. An Experimental study of the
%   collapse of liquid columns on a rigid horizontal plane. Phil. Trans. R. Soc.
%   Lond. A 1952; 244(882):312–324, 1952, DOI: 10.1098/rsta.1952.0006.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       rho1        scalar {998}           Density of the water
%       rho2        scalar {1.205}         Density of the surrounding air
%       miu1        scalar {0.001}         Viscosity of the water
%       miu2        scalar {1.983e-5}      Viscosity of the surrounding air
%       gy          scalar {9.82}          Gravitational constant
%       a           scalar {0.05715}       Width and height of the water column
%       igrid       scalar 1/{0}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {a/12}          Max grid cell size
%       dt          scalar {}              Time step size
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
  'rho1',     998; ...
  'rho2',     1.205; ...
  'miu1',     0.001; ...
  'miu2',     1.983e-5; ...
  'gy',       9.82; ...
  'a',        0.05715; ...
  'igrid',    0; ...
  'hmax',     0.05715/18; ...
  'dt',       {}; ...
  'sf_u',     'sflag2';  ...
  'sf_p',     'sflag1'; ...
  'sf_c',     'sflag1'; ...
  'iplot',    1; ...
  'itest',    0; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;
if( ~got.dt )
  TSTEP = 0.05/sqrt(opt.gy/opt.a);
else
  TSTEP = opt.dt;
end
TMAX= 3.0/sqrt(opt.gy/opt.a);
MAXNIT = 30;
if( opt.itest )
  TMAX = TSTEP;
  MAXNIT = 1;
end

% Geometry and grid generation.
l = 6*opt.a;
h = 1.5*opt.a;
fea.grid = rectgrid(round(l/opt.hmax),round(h/opt.hmax),[0 l;0 h]);
if( opt.igrid==1 )
  fea.grid = quad2tri( fea.grid );
end
fea.sdim = { 'x' 'y' };   % Coordinate names.
n_bdr = max(fea.grid.b(3,:)) + 1;   % Increment number of boundaries.
fea.grid.b(3,1) = n_bdr;


% Problem definition.
dh    = num2str(1.5*opt.hmax);      % Smoothing region width.
dw    = ['(c/',dh,'/(sqrt(cx^2+cy^2+eps)))'];
smhs  = ['((0.5*(1+',dw,'+1/pi*sin(pi*',dw,')))*(',dw,'>-1)*(',dw,'<1)+(',dw,'>=1))'];   % Smooth heaviside function.
rho   = [num2str(opt.rho1),'+',num2str(opt.rho2-opt.rho1),'*',smhs];
miu   = [num2str(opt.miu1),'+',num2str(opt.miu2-opt.miu1),'*',smhs,'+2e2*',num2str(opt.hmax),'*sqrt(u^2+v^2+1e-2)'];

% Add Navier-Stokes equations physics mode.
fea   = addphys(fea,@navierstokes);
fea.phys.ns.eqn.coef{1,end} = { rho };
fea.phys.ns.eqn.coef{2,end} = { miu };
fea.phys.ns.eqn.coef{4,end} = { ['-(',rho,')*',num2str(opt.gy)] };
fea.phys.ns.bdr.sel(n_bdr)  = 4;   % Set pressure to zero on last boundary segment.
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };

fea = addphys(fea,@convectiondiffusion);
fea.phys.cd.sfun            = { opt.sf_c };
fea.phys.cd.eqn.coef{2,4}   = { 0.001 };   % Set diffusion coefficient.
fea.phys.cd.eqn.coef{3,4}   = { 'u' };    % Convection velocity in x-direction.
fea.phys.cd.eqn.coef{4,4}   = { 'v' };    % Convection velocity in y-direction.


% Parse physics modes.
fea = parsephys(fea);

% Implement slip boundary conditions on left and bottom walls.
fea.bdr.d{1}{1} = [];
fea.bdr.d{1}{5} = [];
fea.bdr.d{2}{4} = [];
fea.bdr.n{1}{1} = 0;
fea.bdr.n{1}{5} = 0;
fea.bdr.n{2}{5} = 0;

% Parse problem.
fea = parseprob(fea);


% Call to time-dependent solver.
init = { '0', '0', '0', ['(x-',num2str(opt.a),')*((x-',num2str(opt.a),')>(y-',num2str(opt.a),'))+(y-',num2str(opt.a),')*((x-',num2str(opt.a),')<=(y-',num2str(opt.a),'))'] };
fea.sol.u = solvetime( fea, ...
                       'fid',     fid, ...
                       'tmax',    TMAX, ...
                       'tstep',   TSTEP, ...
                       'icub',    3, ...
                       'maxnit',  MAXNIT, ...
                       'init',    init, ...
                       'ischeme', 3 );


% Postprocessing.
if ( opt.iplot>0 )
  figure('position',[0 0 1920 1080])
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
  postplot(fea,'surfexpr',rho)
  % postplot(fea,'isoexpr','c','isolev',[0 0],'arrowexpr',{'u','v'})
  plot([0 l l 0 0],[0 0 h h 0],'k')
  title('Interface')
end


out.err  = [];
out.pass = [];
if ( nargout==0 )
  clear fea out
end
