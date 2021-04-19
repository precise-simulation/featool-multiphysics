function [ fea, out ] = ex_electrostatics2( varargin )
%EX_ELECTROSTATICS2 Axisymmetric model of a spherical capacitor.
%
%   [ FEA, OUT ] = EX_ELECTROSTATICS2( VARARGIN ) Axisymmetric model of a spherical capacitor
%   compared with an analytical solution.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       sfun        string {sflag2}        Shape function for pressure
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'iplot',    1;
            'tol',      5e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Model constants and parameters.
r1   = 0.003;
r2   = 0.01;
r3   = 0.012;
eps0 = 8.85e-12;
epsr = 3.9;
q0   = 6e-11;


% Geometry and grid generation.
fea.sdim = { 'r' 'z' };
fea.geom.objects = { gobj_circle([0 0],r1,'C1') ...
                     gobj_circle([0 0],r2,'C2') ...
                     gobj_circle([0 0],r3,'C3') ...
                     gobj_rectangle(-1.1*r1,0,-1.1*r3,1.1*r3,'R1') ...
                     gobj_rectangle(-1.1*r2,0,-1.1*r3,1.1*r3,'R2') ...
                     gobj_rectangle(-1.1*r3,0,-1.1*r3,1.1*r3,'R3') };
fea = geom_apply_formula( fea, 'C1-R1' );
fea = geom_apply_formula( fea, 'C2-R2' );
fea = geom_apply_formula( fea, 'C3-R3' );

geom_fix = false;
if( geom_fix )   % Manual geometry fix.
  fea.geom.objects = [ fea.geom.objects, fea.geom.objects(2) ];
  fea.geom.objects = [ fea.geom.objects, fea.geom.objects(1) ];
  fea.geom.objects{4}.tag = 'CS4';
  fea.geom.objects{5}.tag = 'CS5';
  fea = geom_apply_formula( fea, 'CS2-CS1' );
  fea = geom_apply_formula( fea, 'CS3-CS4' );
end

fea.grid = gridgen( fea, 'hmax', (r3-r2)/2, 'fid', opt.fid );


% Problem definition.
fea.const = { 'sigma'  { 6e7     0      6e7   } ;
              'eps0'   { eps0    eps0   eps0  } ;
              'epsr'   { 1       epsr   1     } ;
              'tscale' { 1e-17   1e-17  1e-17 } ;
              'rho'    { q0*3/4/pi/(r1^3) 0 -q0*3/4/pi/(r3^3-r2^3) } };

fea = addphys( fea, {@electrostatics 1} );
fea.phys.es.eqn.coef{1,end} = { 'sigma+epsr*eps0/tscale' };
fea.phys.es.eqn.coef{4,end} = { 'rho/tscale' };
fea.phys.es.sfun            = { opt.sfun };
n_bdr = max(fea.grid.b(3,:));
fea.phys.es.bdr.sel(fea.phys.es.bdr.sel>0) = 4;
fea.phys.es.bdr.sel(findbdr(fea,['sqrt(r.^2+z.^2)>',num2str((r2+r3)/2)])) = 2;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'sigma*sqrt(Vr^2+Vz^2)', 'isoexpr', 'V', 'isocolor', 'w' )
  title( 'Surface: current density, contour: electric potential' )
end


% Error checking.
[Vmin,Vmax] = minmaxsubd( 'V', fea );
eeng = intsubd( 'eps0*epsr*(Vr^2+Vz^2)*pi*r', fea );

cap1 = q0/( Vmax - Vmin );
cap2 = q0^2/(2*eeng);
cap_ref = 4*pi*epsr*eps0/( 1/r1 - 1/r2 );

out.err  = [ abs(cap_ref-cap1)/cap_ref ...
             abs(cap_ref-cap2)/cap_ref ];
out.pass = all(out.err<opt.tol);


if ( nargout==0 )
  clear fea out
end
