function [ fea, out ] = ex_magnetostatics5( varargin )
%EX_MAGNETOSTATICS5 Magnetic field test model.
%
%   [ FEA, OUT ] = EX_MAGNETOSTATICS5( VARARGIN ) Magnetic field test model.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       sfun        string {sflag2}        Shape function for pressure
%       hmax        scalar {0.1}           Grid size
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'sfun',     'sflag2';
            'hmax',     0.1;
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.grid = rectgrid( round(1/opt.hmax) );
fea.grid.s( selcells(fea,'y<=(0.5+sqrt(eps))') ) = 2;


% Problem definition.
fea = addphys( fea, @magnetostatics );
fea.phys.ms.eqn.coef{3,end} = { 1 0 };
fea.phys.ms.eqn.coef{4,end} = { 0 1 };
fea.phys.ms.sfun = { opt.sfun };


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 2, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', fea.phys.ms.eqn.vars{2,2}, ...
            'isoexpr', fea.phys.ms.eqn.vars{2,2}, 'isolev', 25, ...
            'arrowexpr', fea.phys.ms.eqn.vars{8,2}, 'arrowcolor', 'w', 'arrowspacing', [45 30] )
  title( 'Magnetic field' )
end


% Error checking.
Az1  = intsubd( fea.phys.ms.eqn.vars{1,2}, fea, 1 );
Az2  = intsubd( fea.phys.ms.eqn.vars{1,2}, fea, 2 );
Mf1  = intsubd( fea.phys.ms.eqn.vars{2,2}, fea, 1 );
Mf2  = intsubd( fea.phys.ms.eqn.vars{2,2}, fea, 2 );
Scb1 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 1 );
Scb2 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 2 );
Scb3 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 3 );
Scb4 = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 4 );
out.err = [ abs(Az1+3.178289e-8)/3.178289e-8 ;
            abs(Az2+3.178289e-8)/3.178289e-8 ;
            abs(Mf1-0.389419)/0.389419 ;
            abs(Mf2-0.520161)/0.520161 ;
            abs(Scb1-0.16233)/0.16233 ;
            abs(Scb2-0.836552)/0.836552 ;
            abs(Scb3+0.83767)/0.83767 ;
            abs(Scb4+0.16343)/0.16343 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
