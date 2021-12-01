function [ fea, out ] = ex_heattransfer2( varargin )
%EX_HEATTRANSFER2 1D Stationary heat transfer with radiation.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER2( VARARGIN ) NAFEMS T2 benchmark
%   example for heat transfer with radiation [1]. The left end of a 0.1 m
%   rod is held at a temperature of 1000 K while the right end is radiating
%   with an emissivity, em = 0.98, and Stefan-Bolzmann constant,
%   sigma = 5.67e-8 Wm^2/K^4.
%
%              +---------- L=0.1m ----------+ T(0.1)?
%          T=1000K                   q_n = em*sigma*(T_amb^4-T^4)
%
%   The rod is made of iron with density 7850 kg/m^3, heat capacity 460
%   J/kgK, and thermal conductivity 55.563 W/mK. The steady state
%   temperature at the left end is sought when the surrounding ambient
%   temperature is 300 K.
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
%       hmax        scalar {0.02}          Grid cell size
%       sfun        string {sflag1}        Finite element shape function
%       solver      string fenics/{}       Use FEniCS or default solver
%       istat       scalar {1}/0           Use stationary (=1), or time dependent solver
%       iplot       scalar {1}/0           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.02;
            'sfun',     'sflag1';
            'solver',   '';
            'istat',    1;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Grid generation.
L        = 0.1;
nx       = round(L/opt.hmax);
fea.grid = linegrid( nx, 0, L );


% Problem definition.
fea.sdim  = { 'x' };                      % Space coordinate name.
fea = addphys( fea, @heattransfer );      % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };          % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.coef{1,end} = 7850;       % Density
fea.phys.ht.eqn.coef{2,end} =  460;       % Heat capacity.
fea.phys.ht.eqn.coef{3,end} =   55.563;   % Thermal conductivity.
fea.phys.ht.eqn.coef{6,end} = { 1000 };   % Initial temperature.

% Boundary conditions.
fea.phys.ht.bdr.sel = [ 1 4 ];
fea.phys.ht.bdr.coef{1,end} = { 1000 [] };
fea.phys.ht.bdr.coef{4,end}{2}{4} = '0.98*5.67e-8';
fea.phys.ht.bdr.coef{4,end}{2}{5} = 300;


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid, ...
                'tstep', 10, 'tmax', 1000, 'ischeme', 2*(~opt.istat) );
else
  if( opt.istat )
    fea.sol.u = solvestat( fea, 'fid', opt.fid, 'init', {'T0_ht'} );
  else
    [fea.sol.u, tlist] = solvetime( fea, 'fid', opt.fid, 'init', {'T0_ht'}, ...
                                    'tmax', 1000, 'tstep', 10 );
  end
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'T', 'axequal', 0 )
  title('Temperature')
  xlabel('x')
  ylabel('T')
end


% Error checking.
T_sol = evalexpr( 'T', 0.1, fea );
T_ref = 926.97;
out.err  = abs(T_sol-T_ref)/T_ref;
out.pass = out.err<6e-4;


if( nargout==0 )
  clear fea out
end
