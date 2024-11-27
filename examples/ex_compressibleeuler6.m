function [ fea, out ] = ex_compressibleeuler6( varargin )
%EX_COMPRESSIBLEEULER6 2D Compressible inviscid flow past a wedge.
%
%   [ FEA, OUT ] = EX_COMPRESSIBLEEULER6( VARARGIN ) This verification
%   test case studies steady inviscid compressible flow past a wedge
%   at an incident angle of 15 degrees. As the supersonic flow (Ma =
%   2.5) hits the wedge a sharp oblique shock wave is formed,
%   resulting in a reduced downstream flow velocity. The simulation
%   uses the inviscid compressible Euler equations to model the flow,
%   adaptively refining the mesh, and using FEM-TVD upwinding to
%   stabilize the solution and resolve shock discontinuities.
%
%   The angle of the shock wave and downstream Mach number can be
%   determined using oblique shock theory, Ma = 1.873526 at an angle
%   of 36.9449 degrees, and also comparing to results from the
%   NASA/NPARC CFD Verification and Validation Database and from the
%   Wind-US CFD code [1].
%
%   Reference:
%
%   [1] NPARC Alliance, Computational Fluid Dynamics (CFD)
%   Verification and Validation Web Site
%   https://www.grc.nasa.gov/WWW/wind/valid/wedge/wedge.html
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

% Copyright 2013-2024 Precise Simulation, Ltd.


cOptDef = { 'hmax',     0.05;
            'sfun',     'sflag1';
            'solver',   '';
            'endTime',  2;
            'iplot',    1;
            'tol',      0.03;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


fea.sdim = { 'x', 'y' };
fea.geom.objects = { gobj_rectangle( 0, 2, 0, 1, 'R1' ) };
fea = geom_split_object( fea, { 'R1' }, [1 0], [1 atan(pi*15/180)] );
fea = geom_split_object( fea, { 'SP1' }, [1 0], [1 atan(pi*(15+36.9449)/180)] );
fea.geom.objects(3) = [];
hmaxb = opt.hmax * ones(1,7); hmaxb(3) = opt.hmax / 50;
fea.grid = gridgen( fea, 'gridgen', 'default', 'hmax', opt.hmax, 'hmaxb', hmaxb, 'fid', opt.fid );

fea = addphys(fea,@compressibleeuler);

Ma    = 2.5;
rho0  = 1;
p0    = 1;
u0    = Ma*sqrt(1.4*p0/rho0);
v0    = 0;
init0 = {rho0, u0, v0, p0};
[fea.phys.ee.eqn.coef{5,end}{:}] = deal(rho0);
[fea.phys.ee.eqn.coef{6,end}{:}] = deal(u0);
[fea.phys.ee.eqn.coef{7,end}{:}] = deal(v0);
[fea.phys.ee.eqn.coef{8,end}{:}] = deal(p0);

i_in = findbdr( fea, 'x<=sqrt(eps)' );
i_out = findbdr( fea, 'x>=2-sqrt(eps)' );
fea.phys.ee.bdr.sel(i_in) = 1;
fea.phys.ee.bdr.sel(i_out) = 2;
fea.phys.ee.bdr.coef{1,end}{1,i_in} = rho0;
fea.phys.ee.bdr.coef{1,end}{2,i_in} = u0;
fea.phys.ee.bdr.coef{1,end}{3,i_in} = v0;
fea.phys.ee.bdr.coef{1,end}{4,i_in} = p0;

fea.phys.ee.prop.artstab.id = 0;
fea.phys.ee.prop.artstab.sd = 0;
fea.phys.ee.prop.artstab.iupw = 1;

fea = parsephys(fea);
fea = parseprob(fea);


if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'deltaT', 0.005, 'endTime', opt.endTime, 'maxCo', 0.5, 'fid', fid, 'logfid', logfid, 'nproc', 1 );
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
  postplot( fea, 'surfexpr', s_Ma, 'isoexpr', s_Ma )
  title(fea.phys.ee.eqn.vars{end-1,1})
end


% Error checking.
out.err = abs(evalexpr('sqrt(u^2+v^2)/sqrt(ga_ee*p/rho)',[1.9;0.4],fea) - 1.873526) / 1.873526;
out.pass = all(out.err<opt.tol);

if( nargout==0 )
  clear fea out
end
