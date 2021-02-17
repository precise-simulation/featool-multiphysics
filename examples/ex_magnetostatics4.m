function [ fea, out ] = ex_magnetostatics4( varargin )
%EX_MAGNETOSTATICS4 2D Axisymmetric cylindrical magnet example.
%
%   [ FEA, OUT ] = EX_MAGNETOSTATICS4( VARARGIN ) 2D Axisymmetric magnetostatic
%   test example for a cylindrical magnet with no electrical currents.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       sfun        string {sflag2}        Shape function for pressure
%       hmax        scalar {0.0025}        Max grid cell size
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'iplot',    1;
            'hmax',     0.0025;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'r' 'z' };
fea.geom.objects = { gobj_rectangle( 0, 0.03,  -0.03, 0.03, 'R1' ) ...
                     gobj_rectangle( 0, 0.005, -0.01, 0.01, 'R2' ) };

fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', opt.fid );


% Problem definition.
fea = addphys( fea, {@magnetostatics 1} );
fea.phys.ms.eqn.coef{4,end} = { 0 1 };
fea.phys.ms.sfun            = { opt.sfun };


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', fea.phys.ms.eqn.vars{1,2}, 'arrowexpr', fea.phys.ms.eqn.vars{9,2} )
  title( 'Magnetic potential (surface), flux density (arrow)' )
end


% Error checking.
[Atmin,Atmax] = minmaxsubd( fea.phys.ms.eqn.vars{1,2}, fea );
i_rbdr = findbdr( fea, ['r>=0.03-',num2str(sqrt(eps))] );
Scb2 = intbdr( fea.phys.ms.bdr.vars{2,2}, fea, i_rbdr );
out.err = [ abs(Atmax-2.824e-9)/2.814e-9 ;
            abs(Scb2+4.876998e-4)/4.876998e-4 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
