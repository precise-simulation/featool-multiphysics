function [ fea, out ] = ex_linearelasticity5( varargin )
%EX_LINEARELASTICITY5 Vibration of a square plate.
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY5( VARARGIN ) Vibration of a
%   square plate (NAFEMS FV52 Benchmark).
%
%   Reference:
%
%   [1] National Agency for Finite Element Methods and Standards. The
%   Standard NAFEMS Benchmarks. Rev. 3. United Kingdom: NAFEMS,
%   October 1990.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 0/{1}           Cell type (>0=hexahedral, <0=tetrahedral)
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
            'tol',      0.02;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


E   = 200e9;
nu  = 0.3;
rho = 8000;


% Geometry definition.
fea.sdim = {'x' 'y' 'z'};
gobj = gobj_block( 0, 10, 0, 10, -0.5, 0.5, 'B1' );
fea.geom.objects = { gobj };


fea.grid = blockgrid(abs(opt.igrid)*10,abs(opt.igrid)*10,abs(opt.igrid),[0,10;0,10;-0.5,0.5]);
if( opt.igrid<0 )
  fea.grid = hex2tet( fea.grid );
end


% Equations and problem definition.
fea = addphys( fea, @linearelasticity );
fea.phys.el.eqn.coef{1,end} = { nu  };
fea.phys.el.eqn.coef{2,end} = { E   };
fea.phys.el.eqn.coef{3,end} = { rho };
fea.phys.el.sfun            = { opt.sfun, opt.sfun, opt.sfun };


% Set/constrain w = 0 on lower edges (z = -0.5).
edg = [];
[be,e,ev] = gridbdre( fea.grid.b, fea.grid.c );
n_bdre = max(be(end,:));
for i_bdre=1:n_bdre
  ix = be(end,:) == i_bdre;
  ie = be(4,ix);
  iv = unique([ev(ie,1);ev(ie,2)]);
  pz = fea.grid.p(3,iv);
  if( all( pz <= -0.5+sqrt(eps)) )
    edg_i.type  = 'constraint';
    edg_i.index = i_bdre;
    edg_i.dvar  = 3;
    edg_i.expr  = 0;

    edg = [ edg, edg_i ];
  end
end
fea.edg = edg;


% Solve problem.
fea = parsephys( fea );
fea = parseprob( fea );

[fea.sol.u,fea.sol.l] = solveeig( fea, 'neigs', 10, 'fid', fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', 'solnum', 4 )
end


out = [];
f = sqrt(max(0,fea.sol.l))/(2*pi);
f_ref = [0;0;0;45.897;109.44;109.44;167.89;193.59;206.19;206.19];
out.err  = norm(f_ref-f)/norm(f_ref);
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
