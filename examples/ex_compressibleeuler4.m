function [ fea, out ] = ex_compressibleeuler4( varargin )
%EX_COMPRESSIBLEEULER4 2D Compressible inviscid flow over a bump.
%
%   [ FEA, OUT ] = EX_COMPRESSIBLEEULER4( VARARGIN ) Sets up and
%   solves a stationary 2D compressible Euler equation problem for
%   supersonic flow over a cylindrical bump. Reference:
%
%   [1] Lynn J.F., van Leer B., Lee D. (1997) Multigrid solution of
%   the euler equations with local preconditioning. In: Kutler P.,
%   Flores J., Chattot JJ. (eds) Fifteenth International Conference on
%   Numerical Methods in Fluid Dynamics. Lecture Notes in Physics,
%   vol 490. Springer, Berlin, Heidelberg.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.05}          Max grid cell size
%       sfun        string {sflag1}        Shape function
%       solver      string openfoam/su2/{} Use OpenFOAM, SU2, or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.05;
            'sfun',     'sflag1';
            'solver',   '';
            'iplot',    1;
            'tol',      0.1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


fea.sdim = { 'x', 'y' };
r = (0.5^2/0.042+0.042)/2;
fea.geom.objects = { gobj_rectangle(-1,4,0,2), gobj_circle([.5,0.042-r],r) };
fea = geom_apply_formula(fea,'R1-C1');
fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);

fea = addphys(fea,@compressibleeuler);

Ma    = 1.4;
rho0  = 1;
p0    = 1;
u0    = Ma*sqrt(1.4*p0/rho0);
v0    = 0;
init0 = {rho0, u0, v0, p0};
fea.phys.ee.eqn.coef{5,end}{1} = rho0;
fea.phys.ee.eqn.coef{6,end}{1} = u0;
fea.phys.ee.eqn.coef{7,end}{1} = v0;
fea.phys.ee.eqn.coef{8,end}{1} = p0;

i_in = findbdr( fea, 'x<=-1+sqrt(eps)' );
i_out = findbdr( fea, 'x>=4-sqrt(eps)' );
fea.phys.ee.bdr.sel(i_in) = 1;
fea.phys.ee.bdr.sel(i_out) = 2;
fea.phys.ee.bdr.coef{1,end}{1,i_in} = rho0;
fea.phys.ee.bdr.coef{1,end}{2,i_in} = u0;
fea.phys.ee.bdr.coef{1,end}{3,i_in} = v0;
fea.phys.ee.bdr.coef{1,end}{4,i_in} = p0;

fea = parsephys(fea);
fea = parseprob(fea);


if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'deltaT', 0.01, 'endTime', 5, 'fid', fid, 'logfid', logfid );
  fid = logfid;
elseif( strcmp(opt.solver,'su2') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = su2( fea, 'upwind', 'jst', 'fid', fid, 'logfid', logfid, 'nproc', 1 );
  fid = logfid;
else
  fea.sol.u = solvestat( fea, 'init', init0, 'maxnit', 50, 'nlrlx', 0.9, 'fid', fid );
end



% Postprocessing.
s_Ma = fea.phys.ee.eqn.vars{end-1,2};
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', s_Ma )
  title(fea.phys.ee.eqn.vars{end-1,1})
end


% Error checking.
[Ma_min,Ma_max] = minmaxsubd( s_Ma, fea );
out.err = [ abs(Ma_min-1.0)/1.0, ...
            abs(Ma_max-1.8)/1.8 ];
out.pass = all(out.err<opt.tol);

if( nargout==0 )
  clear fea out
end
