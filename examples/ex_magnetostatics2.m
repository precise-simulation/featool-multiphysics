function [ fea, out ] = ex_magnetostatics2( varargin )
%EX_MAGNETOSTATICS2 Magnetic field around a horseshoe magnet.
%
%   [ FEA, OUT ] = EX_MAGNETOSTATICS2( VARARGIN ) Magnetic field around a horseshoe magnet.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       sfun        string {sflag2}        Shape function for pressure
%       hmax        scalar {0.01}          Grid size
%       iorient     scalar 0/{1,2,3}       Magnet orientation (top, right, bottom, left)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'hmax',     0.01;
            'iorient',  0;
            'iplot',    1;
            'tol',      1e-1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.geom.objects = { gobj_circle([0 0],0.05,'C1'), ...
                     gobj_circle([0 0],0.025,'C2'), ...
                     gobj_rectangle(-0.06,0.06,0,0.06,'R1'), ...
                     gobj_rectangle(-0.05,-0.025,0,0.06,'R2'), ...
                     gobj_rectangle(0.025,0.05,0,0.06,'R3'), ...
                     gobj_rectangle(-0.15,0.15,-0.2,0.2,'R4') };
fea = geom_apply_formula( fea, 'C1-C2-R1' );
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );
fea.grid = gridrotate( fea.grid, -pi/2*opt.iorient );


% Problem definition.
fea = addphys( fea, @magnetostatics );
switch( opt.iorient )
  case 0
    fea.phys.ms.eqn.coef{4,end} = { 1 -1 0 0 };
  case 1
    fea.phys.ms.eqn.coef{3,end} = { 1 -1 0 0 };
  case 2
    fea.phys.ms.eqn.coef{4,end} = { -1 1 0 0 };
  case 3
    fea.phys.ms.eqn.coef{3,end} = { -1 1 0 0 };
end
fea.phys.ms.sfun = { opt.sfun };


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', 'Az', ...
            'isoexpr', 'Az', 'isolev', 25, ...
            'arrowexpr', fea.phys.ms.eqn.vars{9,2}, 'arrowcolor', 'w', 'arrowspacing', [45 30] )
  title( 'Magnetic potential (surface, iso), and flux density (arrows) ' )
end


% Error checking.
Az   = intsubd( fea.phys.ms.eqn.vars{1,2}, fea );
Mf   = intsubd( fea.phys.ms.eqn.vars{2,2}, fea );
gAzb = intbdr(  fea.phys.ms.eqn.vars{5,2}, fea, 1:4 );
Scb1 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 1 );
Scb2 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 2 );
Scb3 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 3 );
Scb4 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 4 );
out.err = [ abs(Az+9.637833e-11)/9.637833e-11;
            abs(Mf-0.005108)/0.005108;
            abs(gAzb-1.15859e-8)/1.15859e-8;
            abs(Scb1-0.001326)/0.001326;
            abs(Scb2+0.001865)/0.001865;
            abs(Scb3-0.002403)/0.002403;
            abs(Scb4+0.001865)/0.001865 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
