function [ fea, out ] = ex_acconduction2( varargin )
%EX_ACCONDUCTION2 Dissipation factor in a cylindrical capacitor.
%
%   [ FEA, OUT ] = EX_ACCONDUCTION2( VARARGIN ) Calculate the dissipation factor in a cylindrical capacitor.
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

% Copyright 2013-2025 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'r', 'z' };
fea.geom.objects = { gobj_circle( [0, 0], 25e-3, 'airo' ), ...
                     gobj_circle( [0, 0], 4e-3, 'airi' ), ...
                     gobj_polygon( [1e-3, -1.5e-3; 1.5e-3, -1.5e-3; 1.5e-3, -1.35e-3; 1.5e-3, -1.25e-3; 1.5e-3, 1.5e-3; 1e-3, 1.5e-3; 1e-3, 1.4e-3], 'ceramic' ), ...
                     gobj_rectangle( 0, -25e-3, -25e-3, 25e-3, 'left1' ), ...
                     gobj_rectangle( 0, -25e-3, -25e-3, 25e-3, 'left2' ) };
fea = geom_apply_formula( fea, 'airo - left1' );
fea = geom_apply_formula( fea, 'airi - left2' );

fea.grid = gridgen( fea, 'hmax', [0.05, 2, 0.25]*1e-3, 'fid', opt.fid );


% Problem definition.
fea = addphys( fea, {@conductivemediadc, 1} );
fea.phys.dc.eqn.coef{2,end} = { 's_coef' };
fea.phys.dc.sfun            = { opt.sfun };

fea.phys.dc.bdr.sel([6,7,9,12]) = -3;
fea.phys.dc.bdr.coefi{3,end}{6} = 5;
fea.phys.dc.bdr.coefi{3,end}{7} = 5;
fea.phys.dc.bdr.coefi{3,end}{9} = -5;
fea.phys.dc.bdr.coefi{3,end}{12} = 5;

fea.expr = {'freq', 1e3;
            'epsilon0', 8.854187817e-12;
            'epsilon', '1*epsilon0 6 6';
            'sigma', {1.e-8, 0, 0};
            'omega', '2*pi*freq';
            's_coef', 'epsilon-(i*sigma/omega)';
            'jr', '-sigma*Vr';
            'jz', '-sigma*Vz';
            'jdr', '-i*omega*epsilon*Vr';
            'jdz', '-i*omega*epsilon*Vz'};

% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'abs(V)/sqrt(2)', 'isoexpr', 'abs(V)/sqrt(2)' )
  axis([-.0018, .0043, -.0023, .0024])
  title( 'Electric potential (RMS)' )
end


% Error checking.
P = 2*pi*0.5*intsubd( '-Vr*jr-Vz*jz', fea, 1 ) * 1;  % W/10 mm^3
P_ref = 2.45e-8;
Pr = 2*pi*0.5*intsubd( '-Vr*jdr-Vz*jdz', fea, 1 ) * 1;  % W/10 mm^3
Pr_ref = 81.7e-8;

out.err = abs([ (P - P_ref)/P_ref ]);
out.pass = all( out.err < opt.tol );

if( ~nargout )
  clear fea out
end
