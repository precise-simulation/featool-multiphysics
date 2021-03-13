function [ fea, out ] = ex_axistressstrain5( varargin )
%EX_AXISTRESSSTRAIN5 Axisymmetric vibration modes of a hollow cylinder.
%
%   [ FEA, OUT ] = EX_AXISTRESSSTRAIN5( VARARGIN ) Axisymmetric
%   vibration modes of a hollow cylinder (NAFEMS Free Vibration Benchmark 41).
%
%   Reference:
%
%   [1] F. Abassian, D.J. Dawswell, and N.C. Knowles, Free Vibration
%   Benchmarks, Volume 3, NAFEMS, Glasgow, 1987.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 0/{2}           Cell type (>0=quadrilaterals, <0=triangles)
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'igrid',    2;
            'sfun',     'sflag1';
            'iplot',    1;
            'tol',      0.01;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


E   = 2e11;
nu  = 0.3;
rho = 8000;


% Geometry definition.
fea.sdim = {'r' 'z'};
gobj = gobj_rectangle( 1.8, 2.2, 0, 10, 'R1' );
fea.geom.objects = { gobj };


fea.grid = rectgrid( abs(opt.igrid)*2, abs(opt.igrid)*50, [ 1.8, 2.2; 0, 10 ]);
if( opt.igrid<0 )
  fea.grid = quad2tri( fea.grid );
end


% Equations and problem definition.
fea = addphys( fea, @axistressstrain );
fea.phys.css.eqn.coef{1,end} = { nu  };
fea.phys.css.eqn.coef{2,end} = { E   };
fea.phys.css.eqn.coef{3,end} = { rho };
fea.phys.css.sfun            = { opt.sfun, opt.sfun };


% Solve problem.
fea = parsephys( fea );
fea = parseprob( fea );

[fea.sol.u,fea.sol.l] = solveeig( fea, 'fid', fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'sqrt((r*u)^2+w^2)', 'solnum', 2 )
end


out = [];
f = sqrt(max(0,fea.sol.l))/(2*pi);
f_ref = [ 0; 243.773387; 378.534723; 394.046384; 397.467494; 405.041753 ];
out.err  = norm(f_ref-f)/norm(f_ref);
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
