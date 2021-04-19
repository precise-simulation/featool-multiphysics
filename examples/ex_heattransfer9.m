function [ fea, out ] = ex_heattransfer9( varargin )
%EX_HEATTRANSFER9 1D Transient heat diffusion with a point source.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER9( VARARGIN ) Transient heat
%   diffusion problem with a point source at one end and analytic
%   solution.
%
%              +---------- L=20m ----------+ dT/dn = 0
%            Q(x=0) = 1   T(t=0) = 0
%
%   References:
%
%   [1] Carslaw HS, Jaeger JC.Conduction of heat in solids, 2ndEd.,
%   Oxford at the Clarendon Press, 1959.
%
%   [2] Strang G, Fix G. An analysis of the Finite Element Method, 2nd Ed.,
%   Wellesley-Cambridge Press, 2008.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       nel         scalar {100}           Number of grid cells
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       ischeme     scalar {2}/1/3         Time stepping scheme
%       imass       scalar {2}/1/3/4       Mass matrix lumping
%       tmax        scalar {0.2}           Maximum time
%       tstep       scalar {0.01}          Time step size
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'nel',      100;
            'sfun',     'sflag1';
            'solver',   '';
            'ischeme',  2;
            'imass',    2;
            'tmax',     1;
            'tstep',    0.01;
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry and grid.
L = 20;

fea.sdim = { 'x' };
fea.geom.objects = { gobj_line(0,L) };
fea.grid = linegrid( opt.nel, 0, L );


% Physics definition.
rho = 1;
cp  = 1;
k   = 1;
Q   = 0;
T0  = 0;

fea = addphys( fea, @heattransfer );
fea.phys.ht.eqn.coef{1,end} = { rho };
fea.phys.ht.eqn.coef{2,end} = {  cp };
fea.phys.ht.eqn.coef{3,end} = {   k };
fea.phys.ht.eqn.coef{5,end} = {   Q };
fea.phys.ht.eqn.coef{6,end} = {  T0 };
fea.phys.ht.sfun = { opt.sfun };

% Point source.
fea.pnt.type  = 'source';
fea.pnt.index = 0;
fea.pnt.dvar  = 1;
fea.pnt.expr  = 1;

fea = parsephys( fea );
fea = parseprob( fea );


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', opt.tstep, 'tmax', opt.tmax, 'ischeme', opt.ischeme );
else
  [fea.sol.u,fea.sol.t] = solvetime( fea, ...
                                     'tstep',   opt.tstep, ...
                                     'tmax',    opt.tmax, ...
                                     'ischeme', opt.ischeme, ...
                                     'imass',   opt.imass, ...
                                     'init',    { 'T0_ht' }, ...
                                     'fid', opt.fid );
end

% Analytical solution.
x = fea.grid.p;
t = fea.sol.t;
for i=1:length(t)   % Loop over time
  T_ref = 2*(t(i)/pi)^(1/2)*((exp((-x.^2./(4*t(i))))-0.5.*x.*(sqrt(pi/t(i)).*erfc((x./(2*sqrt(t(i))))))));
end


% Error checking.
T_sol = evalexpr( 'T', x, fea );
out.err = norm( T_ref(:) - T_sol(:) );
out.pass = out.err < opt.tol;


% Postprocessing.
if( opt.iplot )
  figure, hold on
  plot( x, T_sol, 'b-',  'linewidth', 2 )
  plot( x, T_ref, 'r--', 'linewidth', 1.5 )

  grid on
  axis normal
  ax = axis;
  ax(1:2) = [0,2];
  axis( ax )

  title( sprintf('Solution at time %g',fea.sol.t(end)) )
  ylabel('T')
  xlabel('x')
  legend( 'Computed solution', 'Analytical solution' )
end


if( nargout==0 )
  clear fea out
end
