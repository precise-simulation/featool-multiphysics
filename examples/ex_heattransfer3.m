function [ fea, out ] = ex_heattransfer3( varargin )
%EX_HEATTRANSFER3 1D Transient heat conduction.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER3( VARARGIN ) NAFEMS T3 benchmark
%   example for one-dimensional transient heat conduction [1]. A 10 cm
%   thick steel plate is assumed to have one surface exposed to a varying
%   temperature, T = 100*sin(pi*t/40), while the temperature at the other
%   side is fixed, T = 0. This problem can be seen as one dimensional
%   along the axis aligned with the thickness.
%
%            +---T(0.02)?---- L=0.02m ------------+
%       T=100*sin(pi*t/40)                        T=0
%
%   The temperature at x = 0.02 m is sought at time t = 36. The material
%   parameters of the plate are, density 7200 kg/m^3, heat capacity
%   440.5 J/kgK, and thermal conductivity 35 W/mK.
%
%   Reference:
%
%      [1] The Standard NAFEMS Benchmarks,
%          The National Agency for Finite Element Standards, UK, 1990.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.005}         Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.005;
            'sfun',     'sflag1';
            'solver',   '';
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Grid generation.
L        = 0.1;
nx       = round(L/opt.hmax);
fea.grid = linegrid( nx, 0, L );


% Problem definition.
fea.sdim  = { 'x' };                  % Space coordinate name.
fea = addphys( fea, @heattransfer );  % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };      % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{1,end} = 7200;   % Density
fea.phys.ht.eqn.coef{2,end} =  440.5; % Heat capacity.
fea.phys.ht.eqn.coef{3,end} =   35;   % Thermal conductivity.

% Boundary conditions.
fea.phys.ht.bdr.sel = [ 1 1 ];
fea.phys.ht.bdr.coef{1,end} = { '100*sin(pi*t/40)' 0 };


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', 0.05, 'tmax', 32, 'ischeme', 2 );
  tlist = fea.sol.t;
else
  [fea.sol.u, tlist] = solvetime( fea, 'fid', opt.fid, 'tmax', 32, 'tstep', 0.05 );
end

% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'T', 'axequal', 'off' )
  title('Temperature distribution at time t = 32 s')
  xlabel('x')
  ylabel('T')

  subplot(1,2,2)
  for isol=1:numel(tlist)
    T(isol) = evalexpr( 'T', 0.02, fea, isol );
  end
  plot(tlist,T)
  title('Temperature at x = 0.02 m')
  xlabel('time')
  ylabel('T')
end


% Error checking.
T_sol = evalexpr( 'T', 0.02, fea );
T_ref = 36.6;
out.err  = abs(T_sol-T_ref)/T_ref;
out.pass = out.err<1e-2;


if( nargout==0 )
  clear fea out
end
