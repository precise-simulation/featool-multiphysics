function [ fea, out ] = ex_electrostatics1( varargin )
%EX_ELECTROSTATICS1 Electrostatic test example.
%
%   [ FEA, OUT ] = EX_ELECTROSTATICS1( VARARGIN ) Electrostatics test example.
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
            'tol',      6e-3;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry and grid generation.
fea.sdim = { 'x' 'y' };
fea.grid = rectgrid(10);


% Problem definition.
fea = addphys( fea, @electrostatics );
fea.phys.es.eqn.coef{1,end} = { 2 };
fea.phys.es.eqn.coef{2,end} = { 3 };
fea.phys.es.eqn.coef{3,end} = { 4 };
fea.phys.es.eqn.coef{4,end} = { 5 };
fea.phys.es.sfun            = { opt.sfun };

fea.phys.es.bdr.sel = [2 3 3 4];
fea.phys.es.bdr.coef{3,end}{1,2} =  1;
fea.phys.es.bdr.coef{3,end}{1,3} = -3;

% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', opt.fid );   % Call to stationary solver.


% Postprocessing.
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', 'V' )
  title( 'Electric potential' )
end


% Error checking.
V = intsubd( 'V', fea );
D = intbdr( fea.phys.es.bdr.vars{2,2}, fea, 1:4 );
out.err = [ abs(V-1.25)/1.25;
            abs(D+5)/5  ];
out.pass = all( out.err < opt.tol );


if ( nargout==0 )
  clear fea out
end
