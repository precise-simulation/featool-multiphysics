function [ fea, out ] = ex_axistressstrain4( varargin )
%EX_AXISTRESSSTRAIN4 Axisymmetric heat induced stress for a brake disk.
%
%   [ FEA, OUT ] = EX_AXISTRESSSTRAIN4( VARARGIN ) Quasi-static axisymmetric
%   simulation of a brake disk under braking. The braking process induces heat through
%   friction with the brake pad, which in turn results in stresses in the brake disk.
%
%   Ref. A. Adamowicz, Axisymmetric FE Model to Analysis of Thermal Stresses in a Brake Disk,
%        Journal of Theoretical and Applied Mechanics, 53, 2, pp. 357-370, Warsaw 2015.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       deltad      scalar {5.5e-3}        1/2 Disk thickness
%       rd          scalar {66e-3}         Disk inner radius
%       rp          scalar {75.5e-3}       Pad inner radius
%       Rd          scalar {113.5e-3}      Disk outer radius
%       sfun        string {sflag2}        Shape function for displacements
%       nlev        scalar {2}             Uniform grid refinement level
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'deltad',   5.5e-3;
            'rd',       66e-3;
            'rp',       75.5e-3;
            'Rd',       113.5e-3;
            'sfun',     'sflag2';
            'nlev',     2;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Brake disk material parameters.
kd     = 1.44e-5;       % Thermal diffsivity [m2/s]
Kd     = 51;            % Thermal conductivity [W/m/K]
rhod   = 7100;          % Density [kg/m3]
cpd    = Kd/kd/rhod;    % Specific heat [J/Kg/K]
E      = 99.97e9;       % Modulus of elasticity [Pa]
nu     = 0.29;          % Poisson's ratio
alpha  = 1.08e-5;       % Coefficient of thermal expansion [1/K]

T0     = 20 + 273.15;   % Initial temperature [K]

% Brake pad material parameters.
kp     = 1.46e-5;       % Thermal diffsivity [m2/s]
Kp     = 34.3;          % Thermal conductivity [W/m/K]
rhop   = 4700;          % Density [kg/m3]
cpp    = Kp/kp/rhop;    % Specific heat [J/Kg/K]

% Simulation parameters.
tmax   = 3.96;          % Maximum time.
nsteps = 200;           % Number of time steps.


% Create computational geometry and grid.
fea.sdim = { 'r' 'z' };
fea.grid = rectgrid( 2^(opt.nlev-1)*45, 2^(opt.nlev-1)*5, [opt.rd opt.Rd;0 opt.deltad] );
[~,ix] = findbdr( fea, ['(r>=',num2str(opt.rp),') & (z>=',num2str(opt.deltad-sqrt(eps)),')'], 0 );
fea.grid.b(3,ix) = max(fea.grid.b(3,:)) + 1;   % Add boundary for brake pad.
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Equations and problem definition.
fea = addphys( fea, @axistressstrain );
fea.phys.css.eqn.coef{1,end} = { nu    };
fea.phys.css.eqn.coef{2,end} = { E     };
fea.phys.css.eqn.coef{6,end} = { alpha };
fea.phys.css.eqn.coef{7,end} = { ['T-',num2str(T0)] };
fea.phys.css.sfun            = { opt.sfun opt.sfun };

fea = addphys( fea, {@heattransfer, 1} );
fea.phys.ht.eqn.coef{1,end}  = { rhod };
fea.phys.ht.eqn.coef{2,end}  = { cpd  };
fea.phys.ht.eqn.coef{3,end}  = { Kd   };
fea.phys.ht.sfun             = { opt.sfun };


% Equation coefficients, constants, and expressions.
fea.expr = { 'eta'    [] [] {} ;
             'f'      [] [] {0.5} ;
             'omega'  [] [] {[]} ;
             'p0'     [] [] {1.47e6} };


% Boundary conditions.
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bctype{2,1} = 1;
fea.phys.css.bdr.coef{1,5} = bctype;


eta   = sqrt(Kd*rhod*cpd) / (sqrt(Kd*rhod*cpd) + sqrt(Kp*rhop*cpp));
f     = 0.5;
p0    = 1.47e6;
qd    = [num2str(eta*f*88.464*p0/(2*pi-0.8)),'*r*(1-t/',num2str(tmax),')'];

fea.phys.ht.bdr.sel(5) = 4;
fea.phys.ht.bdr.coef{4,end}{5}{1} = qd;


% Solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
[fea.sol.u,tlist] = solvetime( fea, 'ischeme', 2, 'icub', 3, 'init', { 0 0 T0 }, 'tmax', tmax, 'tstep', tmax/nsteps, 'fid', fid );


% Postprocessing.
sr  = fea.phys.css.eqn.vars{5,end};
st  = fea.phys.css.eqn.vars{6,end};
svm = fea.phys.css.eqn.vars{1,end};
if( opt.iplot>0 )
  figure
  solnum = numel(tlist);
  niso   = 25;
  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'T-273.15', 'isoexpr', 'T-273.15', 'isolev', niso, 'solnum', solnum )
  title(['Temperature [deg C], t = ',num2str(tlist(solnum)),' s'])
  subplot(2,2,2)
  postplot( fea, 'surfexpr', ['(',sr,')*1e-6'],  'isoexpr', ['(',sr,')*1e-6'],  'isolev', niso, 'solnum', solnum )
  title('radial stress [MPa]')
  subplot(2,2,3)
  postplot( fea, 'surfexpr', ['(',st,')*1e-6'],  'isoexpr', ['(',st,')*1e-6'],  'isolev', niso, 'solnum', solnum )
  title('tangential stress [MPa]')
  subplot(2,2,4)
  postplot( fea, 'surfexpr', ['(',svm,')*1e-6'], 'isoexpr', ['(',svm,')*1e-6'], 'isolev', niso, 'solnum', solnum )
  title('von Mieses stress [MPa]')

  figure
  u_stored = fea.sol.u;
  np = 100;
  times = [0.1 0.2 1 2 3 tmax];
  rz = [ linspace(opt.rd,opt.Rd,np); opt.deltad*ones(1,np) ];
  for i=1:numel(times)
    t_i = times(i);
    ix  = max( find( tlist<t_i ) );
    s   = ( tlist(ix+1) - t_i )/( tlist(ix+1) - tlist(ix) );
    fea.sol.u = u_stored(:,ix) + s*( u_stored(:,ix+1) - u_stored(:,ix) );

    subplot(2,2,1)
    T_i = evalexpr( 'T-273.15', rz, fea );
    plot( rz(1,:), T_i ), hold on, grid on
    text( rz(1,end), T_i(end), [num2str(t_i),' s'] )
    title('Temperature [deg C]')
    axis tight

    subplot(2,2,2)
    sr_i = evalexpr( ['(',sr,')*1e-6'], rz, fea );
    plot( rz(1,:), sr_i ), hold on, grid on
    text( rz(1,end-floor(np/4)), sr_i(end-floor(np/4)), [num2str(t_i),' s'] )
    title('radial stress [MPa]')
    axis tight

    subplot(2,2,3)
    st_i = evalexpr( ['(',st,')*1e-6'], rz, fea );
    plot( rz(1,:), st_i ), hold on, grid on
    text( rz(1,end-floor(np/4)), st_i(end-floor(np/4)), [num2str(t_i),' s'] )
    title('tangential stress [MPa]')
    axis tight

    subplot(2,2,4)
    svm_i = evalexpr( ['(',svm,')*1e-6'], rz, fea );
    plot( rz(1,:), svm_i ), hold on, grid on
    text( rz(1,end-floor(np/4)), svm_i(end-floor(np/4)), [num2str(t_i),' s'] )
    title('von Mieses stress [MPa]')
    axis tight
  end
  fea.sol.u = u_stored;

  figure
  hold on

  npth = 72;
  th   = linspace( 0, 2*pi, npth );

  x = opt.rd*cos(th);
  z = opt.rd*sin(th);
  y = opt.deltad*ones(size(x));
  plot3(x, y,z,'k-')
  plot3(x,-y,z,'k-')

  x = opt.Rd*cos(th);
  z = opt.Rd*sin(th);
  plot3(x, y,z,'k-')
  plot3(x,-y,z,'k-')

  npr = 25;
  g = ringgrid( npr-1, npth-1, opt.rd, opt.Rd );
  r = sqrt( g.p(1,:).^2 + g.p(2,:).^2 );
  T = evalexpr( 'T-273.15', [r;zeros(size(r))], fea );

  h = patch( 'faces', g.c', 'vertices', [g.p(1,:)' zeros(size(g.p,2),1) g.p(2,:)'], 'facevertexcdata', T, 'facecolor', 'interp', 'linestyle', 'none' );
  postplot( fea, 'surfexpr', 'T-273.15', 'colorbar', 'off' )
  fea.grid.p(2,:) = -fea.grid.p(2,:);
  postplot( fea, 'surfexpr', 'T-273.15', 'colorbar', 'off' )
  fea.grid.p(2,:) = -fea.grid.p(2,:);
  view(3)
  axis off
  axis tight
end


out = [];
[Tmin,Tmax] = minmaxsubd( 'T-273.15', fea );
Terr = [ abs((34.179-Tmin)/34.179) abs((88.224-Tmax)/88.224) ];

[stmin,stmax] = minmaxsubd( ['(',st,')*1e-6'], fea );
sterr = [ abs((-21.75-stmin)/-21.75) abs((45.1-stmax)/45.1) ];

[svmmin,svmmax] = minmaxsubd( ['(',svm,')*1e-6'], fea );
svmerr = [ abs((4.478-svmmin)/4.478) abs((46.4-svmmax)/46.4) ];

out.T      = [Tmin Tmax];
out.Terr   = Terr;
out.st     = [stmin stmax];
out.sterr  = sterr;
out.sv     = [svmmin svmmax];
out.svmerr = svmerr;
out.pass   = ~any( [ Terr>0.12 sterr>0.03 svmerr>[0.3 0.01] ] );

if( nargout==0 )
  clear fea out
end
