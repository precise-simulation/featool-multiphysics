function [ fea, out ] = ex_acconduction1( varargin )
%EX_ACCONDUCTION1 Calculate dielectric loss in a plane capacitor.
%
%   [ FEA, OUT ] = EX_ACCONDUCTION1( VARARGIN ) Calculate dielectric
%   loss in a 10 mm^2 plane capacitor with 1MHz AC current.
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

% Copyright 2013-2026 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x', 'y' };
fea.geom.objects = { gobj_rectangle( 0, 0.1e-3, 0, 10e-3 ) };
fea.grid = rectgrid( 2, 200, [0, 0.1e-3; 0, 10e-3] );
% fea.grid = quad2tri(fea.grid);


% Problem definition.
fea = addphys( fea, @conductivemediadc );
fea.phys.dc.eqn.coef{2,end} = { 's_coef' };
fea.phys.dc.sfun            = { opt.sfun };

fea.phys.dc.bdr.sel = [2, 1, 2, 1];
fea.phys.dc.bdr.coef{1,end}{2} = 5;

fea.expr = {'freq', 1e6;
            'epsilon0', 8.854187817e-12;
            'epsilon', '10*epsilon0';
            'sigma', 0.00000556;
            'omega', '2*pi*freq';
            's_coef', 'epsilon-(i*sigma/omega)';
            'jx', '-sigma*Vx';
            'jy', '-sigma*Vy';
            'jdx', '-i*omega*epsilon*Vx';
            'jdy', '-i*omega*epsilon*Vy'};

% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'V', 'arrowexpr', {'jx', 'jy'} )
  axis([0, 0.1e-3, 0, 0.1e-3])
  title( 'Electric potential' )
end


% Error checking.
zlen = 1e-2;  % 10 mm
Ia = intbdr( 'nx*jx+ny*jy', fea, 2 );  % A/m^2 (RMS = /sqrt(2))
Ia_ref = -2.78e-5 / zlen;
Ir = imag(intbdr( 'nx*jdx+ny*jdy', fea, 2 ));
Ir_ref = -0.00278163 / zlen;
P = 0.5*intsubd( '-Vx*jx-Vy*jy', fea ) * zlen;  % W/10 mm^3
P_ref = 6.95002e-5;

out.err = abs([ (Ia - Ia_ref)/Ia_ref, (Ir - Ir_ref)/Ir_ref, (P - P_ref)/P_ref ]);
out.pass = all( out.err < opt.tol );

if( ~nargout )
  clear fea out
end
