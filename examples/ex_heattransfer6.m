function [ fea, out ] = ex_heattransfer6( varargin )
%EX_HEATTRANSFER6 2D axisymmetric heat conduction.
%
%   [ FEA, OUT ] = EX_HEATTRANSFER6( VARARGIN ) NAFEMS benchmark example
%   for heating of a solid cylider with an internal hole.
%
%         _              T=T_amb
%         ^            +---------+
%         |    q_n=0   |         |
%         |            :         |
%       0.14m  q_n=5e5 |         | T=T_amb
%         |            :         |
%         |    q_n=0   |         |
%         v            +---------+
%                        T=T_amb
%                      r=0.02
%                      |<-0.08m->|
%
%   The geometry can be considered axisymmetric and the solid has a thermal conductivity
%   of 52 W/mC, the middle part of the inside of the cylider is heated by 5e5 W/mK. The
%   steady temperature at the point (0.04,0.04) is sought when the surrounding ambient
%   temperature is T_amb = 0 C.
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
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry definition.
gobj = gobj_polygon( [ 0.02 0.1 0.1  0.02 0.02 0.02 ;
                       0    0   0.14 0.14 0.1  0.04 ]' );
fea.geom.objects = { gobj };


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );


% Problem definition.
fea.sdim  = { 'r', 'z' };             % Space coordinate name.
fea = addphys( fea, @heattransfer );  % Add heat transfer physics mode.
fea.phys.ht.sfun = { opt.sfun };      % Set shape function.

% Equation coefficients.
fea.phys.ht.eqn.seqn = '- r*k_ht*(Tr_r + Tz_z) = 0';

fea.phys.ht.eqn.coef{3,end} =    52;       % Thermal conductivity.
fea.phys.ht.eqn.coef{7,end} = { 273.15 };  % Initial temperature.

% Boundary conditions.
fea.phys.ht.bdr.sel = [1 1 1 3 4 3];
fea.phys.ht.bdr.coef{1,end}   = { 273.15 273.15 273.15 [] [] [] };
fea.phys.ht.bdr.coef{4,end}{5}{1} = 'r*5e5';


% Parse physics modes and problem struct.
fea = parsephys(fea);
fea = parseprob(fea);


% Compute solution.
if( strcmp(opt.solver,'fenics') )
  fea = fenics( fea, 'fid', opt.fid );
else
  fea.sol.u = solvestat( fea, 'fid', opt.fid, 'init', {'T0_ht'} );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'T', 'isoexpr', 'T' )
  title('Temperature, T')
end


% Error checking.
T_sol = evalexpr( 'T', [0.04;0.04], fea );
T_ref = 332.97;
out.err  = abs(T_sol-T_ref)/T_ref;
out.pass = out.err<opt.tol;


if( nargout==0 )
  clear fea out
end
