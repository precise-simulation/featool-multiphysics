function [ fea, out ] = ex_linearelasticity6( varargin )
%EX_LINEARELASTICITY6 NAFEMS LE6 skew plate benchmark
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY6( VARARGIN ) Skew plate under
%   normal pressure (NAFEMS LE6 Benchmark).
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
%       hmax        string {0.05}          Grid size
%       sfun        string {sflag2}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',    0.05;
            'sfun',     'sflag2';
            'iplot',    1;
            'tol',      0.1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


E    = 210e9;
nu   = 0.3;
rho  = 7800;
t    = 0.01;
load = -0.7e3;


% Grid definition.
fea.sdim = {'x', 'y', 'z'};
a = cos(pi*30/180);
p = [0 1 1+a a;0 0 0.5 0.5]';
geom_2d.objects{1} = gobj_polygon(p);
grid_2d = gridgen( geom_2d, 'hmax', opt.hmax, 'fid', fid );
fea.grid = gridextrude( grid_2d, max(3,ceil(t/opt.hmax)), t );


% Equations and problem definition.
fea = addphys( fea, @linearelasticity );
fea.phys.el.eqn.coef{1,end} = { nu   };
fea.phys.el.eqn.coef{2,end} = { E    };
fea.phys.el.eqn.coef{3,end} = { rho  };
fea.phys.el.sfun            = { opt.sfun, opt.sfun, opt.sfun };


% Set/constrain w = 0 on sides (boundaries 1-4).
n_bdr = 6;
bctype = mat2cell( zeros(3,n_bdr), [1 1 1], ones(1,n_bdr) );
[bctype{3,[1:4]}] = deal(1);
fea.phys.el.bdr.coef{1,5} = bctype;

% Set/constrain u, v = 0 on edge 2, x = y = 0.
edg(1).index = 2;
edg(1).type = 'constraint';
edg(1).dvar = 1;
edg(1).expr = 0;

edg(2).index = 2;
edg(2).type = 'constraint';
edg(2).dvar = 2;
edg(2).expr = 0;

% Set/constrain v = 0 on edge 1, x = 0, y = 1.
edg(3).index = 1;
edg(3).type = 'constraint';
edg(3).dvar = 2;
edg(3).expr = 0;
fea.edg = edg;

% Apply vertical load to top boundary.
bccoef = mat2cell( zeros(3,n_bdr), [1 1 1], ones(1,n_bdr) );
bccoef{3,6} = load;
fea.phys.el.bdr.coef{1,end} = bccoef;


% Solve problem.
fea = parsephys( fea );
fea = parseprob( fea );

fea.sol.u = solvestat( fea, 'fid', fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', ...
            'deformexpr', {'u', 'v', 'w'} )
end


% Error checking.
out = [];
p_E = [0.933012701892219; 0.25; 0];
ps1_E = evalexpr( fea.phys.el.eqn.vars{12,2}, p_E, fea );
ps2_E = evalexpr( fea.phys.el.eqn.vars{13,2}, p_E, fea );
ps3_E = evalexpr( fea.phys.el.eqn.vars{14,2}, p_E, fea );
ps_E_max = max([ps1_E,ps2_E,ps3_E]);
ps_E_ref = 0.802e6;
out.w_E = evalexpr( 'w', p_E, fea );
out.ps_E = [ps1_E, ps2_E, ps3_E];
out.err  = abs(ps_E_max-ps_E_ref)/ps_E_ref;
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
