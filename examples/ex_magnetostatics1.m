function [ fea, out ] = ex_magnetostatics1( varargin )
%EX_MAGNETOSTATICS1 Magnetostatic test example.
%
%   [ FEA, OUT ] = EX_MAGNETOSTATICS1( VARARGIN ) Magnetostatics test example.
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
            'tol',      1.5e-1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.grid = rectgrid(10);


% Problem definition.
fea = addphys( fea, @magnetostatics );
fea.phys.ms.eqn.coef{1,end} = {  1 };
fea.phys.ms.eqn.coef{2,end} = {  2 };
fea.phys.ms.eqn.coef{3,end} = {  3 };
fea.phys.ms.eqn.coef{4,end} = { -4 };
fea.phys.ms.sfun            = { opt.sfun };

fea.phys.ms.bdr.sel = [2 4 3 1];
fea.phys.ms.bdr.coef{3,end}{1,3} = 6;
fea.phys.ms.bdr.coef{1,end}{1,4} = 'y';


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'Az', 'arrowexpr', fea.phys.ms.eqn.vars{8,2} )
  title( 'Magnetic potential' )
  subplot(1,2,2)
  postplot( fea, 'surfexpr', fea.phys.ms.eqn.vars{5,2}, 'arrowexpr', fea.phys.ms.eqn.vars{9,2} )
  title( 'Magnetic field' )
end


% Error checking.
Az  = intsubd( 'Az', fea );
Mf  = intsubd( fea.phys.ms.eqn.vars{2,2}, fea );
Scb = intbdr(  fea.phys.ms.bdr.vars{2,2}, fea, 1:4 );
out.err = [ abs(Az-3.21)/3.21  ;
            abs(Mf-4.9)/4.9    ;
            abs(Scb+2.02)/2.02 ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
