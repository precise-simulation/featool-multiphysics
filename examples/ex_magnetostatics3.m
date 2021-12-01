function [ fea, out ] = ex_magnetostatics3( varargin )
%EX_MAGNETOSTATICS3 3D Cylindrical magnet example.
%
%   [ FEA, OUT ] = EX_MAGNETOSTATICS3( VARARGIN ) 3D Magnetostatic test example for
%   a cylindrical magnet with no electrical currents.
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
            'tol',      1e-1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' 'z' };
fea.geom.objects = { gobj_cylinder([0 0 -0.03],0.03,0.06,3,'C1') ...
                     gobj_cylinder([0 0 -0.01],0.005,0.02,3,'C2') };

nref = 1;
g1 = cylgrid(nref*3,nref*2,nref*12,0.005,0.06,[0 0 -0.03]',3);
g2 = gridextrude( ringgrid(nref*4,nref*12,0.005,0.03,[0;0],pi/4), nref*12, 0.06 );
g2.p(3,:) = g2.p(3,:) - 0.03;
fea.grid  = gridmerge( g1, 1:4, g2, 1:4 );
fea.grid.s(:) = 1;
fea.grid.s( selcells( fea, '(z>=(-0.01-sqrt(eps))).*(z<=(0.01+sqrt(eps))).*(sqrt(x^2+y^2)<=(0.005+sqrt(eps)))' ) ) = 2;


% Problem definition.
fea = addphys( fea, @magnetostatics );
fea.phys.ms.eqn.coef{4,end} = { 0 1 };
fea.phys.ms.sfun            = { opt.sfun };


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'sliceexpr', fea.dvar{1}, 'arrowexpr', fea.phys.ms.eqn.vars{11,2} )
  title( 'Magnetic potential (slice), flux density (arrow)' )
end


% Error checking.
[Vmmin Vmmax]  = minmaxsubd( fea.dvar{1}, fea );
Mf1 = intsubd( fea.phys.ms.eqn.vars{2,2}, fea, 1 );
Mf2 = intsubd( fea.phys.ms.eqn.vars{2,2}, fea, 2 );
out.err = [ abs(Vmmin+0.002158)/0.002158 ;
            abs(Vmmax-0.002158)/0.002158 ;
            abs(Mf1-2.693024e-6)/2.693024e-6 ;
            abs(Mf2-3.084971e-7)/3.084971e-7 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
