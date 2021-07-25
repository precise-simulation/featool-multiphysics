function [ fea, out ] = ex_heattransfer4( varargin )
%EX_HEATTRANSFER4 2D Heat transfer with convective cooling.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER4( VARARGIN ) NAFEMS T4 benchmark
%   example for two dimensional heat transfer with convective heat flux
%   boundary conditions.
%
%        _            q_n=h*(T_amb-T)
%        ^        +--------+
%        |        |        |
%        |  q_n=0 |        | q_n=h*(T_amb-T)
%       1m        |        |
%        |        |   T(0.6,0.2)?
%        |        |        |
%        v        +--------+
%                    T=100
%
%                 |<-0.6m->|
%
%   A 0.6 by 1 m iron plate, with density 7850 kg/m^3, heat capacity 460
%   J/kgC, and thermal conductivity 52 W/mC, is prescribed a fixed
%   temperature of T = 100 C at the bottom edge. The left side is
%   insulated, and the right and top boundaries exposed to convective
%   cooling with a heat transfer coefficient h = 750 W/m^2K. The steady
%   state temperature at the point (0.6,0.2) is sought when the surrounding
%   ambient temperature is T_amb = 0 C.
%
%   Reference:
%
%      [1] Cameron AD, Casey JA, Simpson GB. Benchmark Tests for Thermal Analysis,
%          The National Agency for Finite Element Standards, UK, 1986.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.025}         Grid cell size
%       igrid       scalar {0}/1/2         Cell type (0=quadrilaterals, 1=triangles,
%                                          2=triangles converted from quadrilaterals)
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       istat       scalar {1}/0           Use stationary (=1), or time dependent solver
%       iplot       scalar {1}/1           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.025;
            'igrid',    0;
            'sfun',     'sflag1';
            'solver',   '';
            'istat',    1;
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry definition.
gobj = gobj_rectangle( 0, 0.6, 0, 1 );
fea.geom.objects = { gobj };


% Grid generation.
switch opt.igrid
  case 0
    fea.grid = rectgrid( round(0.6/opt.hmax), round(1/opt.hmax), [0 0.6;0 1] );
  case 1
    fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );
  case 2
    fea.grid = rectgrid( round(0.6/opt.hmax), round(1/opt.hmax), [0 0.6;0 1] );
    fea.grid = quad2tri( fea.grid, 1 );
end


% Problem definition.
fea.sdim  = { 'x', 'y' };             % Space coordinate name.
fea = addphys( fea, @heattransfer );  % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };      % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{1,end} = 7850;   % Density
fea.phys.ht.eqn.coef{2,end} =  460;   % Heat capacity.
fea.phys.ht.eqn.coef{3,end} =   52;   % Thermal conductivity.
fea.phys.ht.eqn.coef{7,end} = { 0 };  % Initial temperature.

% Boundary conditions.
fea.phys.ht.bdr.sel = [1 4 4 3];
fea.phys.ht.bdr.coef{1,end}   = { 100 [] [] [] };
fea.phys.ht.bdr.coef{4,end}{2}{2} = 750;
fea.phys.ht.bdr.coef{4,end}{3}{2} = 750;


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', 100, 'tmax', 20000, 'ischeme', 2*(~opt.istat) );
else
  if( opt.istat )
    fea.sol.u = solvestat( fea, 'fid', opt.fid, 'init', {'T0_ht'} );
  else
    [fea.sol.u, tlist] = solvetime( fea, 'fid', opt.fid, 'init', {'T0_ht'}, ...
                                    'tmax', 20000, 'tstep', 100, 'toldef', 1e-4, 'maxnit', 5 );
  end
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'T', 'isoexpr', 'T' )
  title('Temperature, T')
end


% Error checking.
T_sol = evalexpr( 'T', [0.6;0.2], fea );
T_ref = 18.3;
out.err  = abs(T_sol-T_ref)/T_ref;
out.pass = out.err<opt.tol;


if( nargout==0 )
  clear fea out
end
