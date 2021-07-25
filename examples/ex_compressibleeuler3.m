function [ fea, out ] = ex_compressibleeuler3( varargin )
%EX_COMPRESSIBLEEULER3 2D Steady reflected shock problem.
%
%   [ FEA, OUT ] = EX_COMPRESSIBLEEULER3( VARARGIN ) Sets up and
%   solves a steady reflected compressible Euler shock problem and
%   compares with the analytical solution.
%
%   Ref. H. W. Liepmann, A. Roshko, Elements of Gas Dynamics, Courier Corporation, 2013.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       lev         scalar {2}             Grid refinement level
%       sfun        string {sflag1}        Shape function
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2, or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'lev',      2;
            'sfun',     'sflag1';
            'solver',   '';
            'iplot',    1;
            'tol',      0.15;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


fea.sdim = { 'x', 'y' };
fea.geom.objects = { gobj_rectangle(0,4.1,0,1) };
lev = ceil(opt.lev);
fea.grid = rectgrid(lev*41,lev*10,[0,4.1;0,1]);
fea.grid = quad2tri(fea.grid,1);

ML  = 2.9;
rL  = 1;
uL  = 2.9;
vL  = 0;
pL  = 0.714286;

th1 = 29*pi/180;
MT  = 2.3781;
rT  = 1.7;
uT  = 2.619334;
vT  = -0.5063;
pT  = 1.52819;

th2 = 23*pi/180;
MR  = 1.94253;
rR  = 2.6872838;
uR  = 2.401499;
vR  = 0;
pR  = 2.93398;

s = '%g-((1-y)>x*%g)*(%g-%g)*(x<%g)-(x>=%g)*(y<(x-%g)*%g)*(%g-%g)';
rref = sprintf( s, rT, atan(th1), rT, rL, 1/atan(th1), 1/atan(th1), 1/atan(th1), atan(th2), rT, rR );
uref = sprintf( s, uT, atan(th1), uT, uL, 1/atan(th1), 1/atan(th1), 1/atan(th1), atan(th2), uT, uR );
vref = sprintf( s, vT, atan(th1), vT, vL, 1/atan(th1), 1/atan(th1), 1/atan(th1), atan(th2), vT, vR );
pref = sprintf( s, pT, atan(th1), pT, pL, 1/atan(th1), 1/atan(th1), 1/atan(th1), atan(th2), pT, pR );

fea = addphys(fea,@compressibleeuler);
fea.phys.ee.prop.artstab.id_coef = 2*fea.phys.ee.prop.artstab.id_coef;
fea.phys.ee.prop.artstab.sd_coef = 2*fea.phys.ee.prop.artstab.sd_coef;

init0 = { rL, uL, vL, pL };
fea.phys.ee.eqn.coef{5,end}{1} = rL;
fea.phys.ee.eqn.coef{6,end}{1} = uL;
fea.phys.ee.eqn.coef{7,end}{1} = vL;
fea.phys.ee.eqn.coef{8,end}{1} = pL;

fea.phys.ee.bdr.sel(2)   = 2;
fea.phys.ee.bdr.sel(3:4) = 1;
fea.phys.ee.bdr.coef{1,end}{1,3} = rT;
fea.phys.ee.bdr.coef{1,end}{2,3} = uT;
fea.phys.ee.bdr.coef{1,end}{3,3} = vT;
fea.phys.ee.bdr.coef{1,end}{4,3} = pT;
fea.phys.ee.bdr.coef{1,end}{1,4} = rL;
fea.phys.ee.bdr.coef{1,end}{2,4} = uL;
fea.phys.ee.bdr.coef{1,end}{3,4} = vL;
fea.phys.ee.bdr.coef{1,end}{4,4} = pL;

fea = parsephys(fea);
fea = parseprob(fea);


if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'deltaT', 0.01, 'endTime', 4, 'maxCo', 0.2, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  fea.sol.u = solvestat( fea, 'init', init0, 'fid', fid );
end


% Postprocessing.
s_Ma = fea.phys.ee.eqn.vars{end-1,2};
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', s_Ma )
  title(fea.phys.ee.eqn.vars{end-1,1})
end


% Error checking.
r = evalexprp( fea.dvar{1}, fea );
u = evalexprp( fea.dvar{2}, fea );
v = evalexprp( fea.dvar{3}, fea );
p = evalexprp( fea.dvar{4}, fea );
r_ref = evalexprp( rref, fea );
u_ref = evalexprp( uref, fea );
v_ref = evalexprp( vref, fea );
p_ref = evalexprp( pref, fea );
out.err = [ sum(abs(r_ref-r))/size(fea.grid.p,2), ...
            sum(abs(u_ref-u))/size(fea.grid.p,2), ...
            sum(abs(v_ref-v))/size(fea.grid.p,2), ...
            sum(abs(p_ref-p))/size(fea.grid.p,2) ];
out.pass = all(out.err<opt.tol);


if( nargout==0 )
  clear fea out
end
