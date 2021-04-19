function [ fea, out ] = ex_compressibleeuler2( varargin )
%EX_COMPRESSIBLEEULER2 2D Steady oblique shock wave.
%
%   [ FEA, OUT ] = EX_COMPRESSIBLEEULER2( VARARGIN ) Sets up and
%   solves a steady 2D compressible Euler equation for a Ma=2 10
%   degree oblique shock wave and compares with the analytical
%   solution.
%
%   Ref. H. W. Liepmann, A. Roshko, Elements of Gas Dynamics, Courier Corporation, 2013.
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
            'tol',      0.05;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


fea.sdim = { 'x', 'y' };
fea.geom.objects = { gobj_rectangle() };
fea.grid = rectgrid(ceil(1/opt.hmax));

gamma = 7/5;

al  = pi/18;   % 10 degrees indicence angle.
Min = 2;
Ma  = @(th) sqrt( 2/( sin(th+al)*cos(th+al)*( (gamma+1)*tan(th) - (gamma-1)*tan(th+al) ) ) ) - Min;
th  = fzero(Ma,atan(29.3/180*pi));   % Shock angle.

rin = 1;
uin =  cos(al);
vin = -sin(al);
pin = (sqrt(uin^2+vin^2)/Min)^2*rin/gamma;

Mout = 1/sin(th) * sqrt( (1+(gamma-1)/2*Min^2*sin(th+al)^2) / ...
                         (gamma*Min^2*sin(th+al)^2-(gamma-1)/2) );
rout = (gamma+1)*Min^2*sin(th+al)^2/( (gamma-1)*Min^2*sin(th+al)^2 + 2 );
pout = pin*(1 + 2*gamma/(gamma + 1)*( Min^2*sin(th+al)^2 - 1 ));
uout = Mout*sqrt(gamma*pout/rout);

rref = sprintf( '%g+(y<x*%g)*(%g-%g)', rin, atan(th), rout, rin );
uref = sprintf( '%g+(y<x*%g)*(%g-%g)', uin, atan(th), uout, uin );
vref = sprintf( '%g+(y<x*%g)*(%g-%g)', vin, atan(th), 0,    vin );
pref = sprintf( '%g+(y<x*%g)*(%g-%g)', pin, atan(th), pout, pin );

fea = addphys(fea,@compressibleeuler);
fea.phys.ee.prop.artstab.id_coef = 2*fea.phys.ee.prop.artstab.id_coef;
fea.phys.ee.prop.artstab.sd_coef = 2*fea.phys.ee.prop.artstab.sd_coef;

init0 = { rin, uin, vin, pin };
fea.phys.ee.eqn.coef{5,end}{1} = rin;
fea.phys.ee.eqn.coef{6,end}{1} = uin;
fea.phys.ee.eqn.coef{7,end}{1} = vin;
fea.phys.ee.eqn.coef{8,end}{1} = pin;

fea.phys.ee.bdr.sel(2)   = 2;
fea.phys.ee.bdr.sel(3:4) = 1;
fea.phys.ee.bdr.coef{1,end}{1,3} = rin;
fea.phys.ee.bdr.coef{1,end}{2,3} = uin;
fea.phys.ee.bdr.coef{1,end}{3,3} = vin;
fea.phys.ee.bdr.coef{1,end}{4,3} = pin;
fea.phys.ee.bdr.coef{1,end}{1,4} = rin;
fea.phys.ee.bdr.coef{1,end}{2,4} = uin;
fea.phys.ee.bdr.coef{1,end}{3,4} = vin;
fea.phys.ee.bdr.coef{1,end}{4,4} = pin;

fea = parsephys(fea);
fea = parseprob(fea);


if( strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'deltaT', 0.01, 'endTime', 10, 'fid', fid, 'logfid', logfid );
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
