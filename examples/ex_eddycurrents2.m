function [ fea, out ] = ex_eddycurrents2( varargin )
%EX_EDDYCURRENTS2 3D Eddy currents test example.
%
%   [ FEA, OUT ] = EX_EDDYCURRENTS2( VARARGIN ) 3D Eddy currents test
%   example for vector elements (Nedelec). Accepts the following
%   property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       icase       integer 1/{1}          Predefined test case
%       hmax        scalar {0.05}          Grid cell size
%       sfun        string {sf_simp_N1}    Vector shape function (Nedelec)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'icase',    1;
            'hmax',     0.05;
            'sfun',     'sf_simp_N1';
            'iplot',    1;
            'tol',      0.01;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry and grid generation.
fea.sdim = {'x','y','z'};
fea.geom.objects = { gobj_rectangle() };
fea.grid = hex2tet( blockgrid(ceil(1/opt.hmax)), 2 );


% Problem definition.
fea = addphys( fea, @customeqn, {'E'} );
fea.phys.ce.eqn.seqn = {'Ec_c + E_t = 0'};
fea.phys.ce.sfun = {opt.sfun};
if( opt.icase==2 )   % Set homogenous Neumann BCs.
  [fea.phys.ce.bdr.coef{5}{:}] = deal(0);
end

% Parse problem.
fea = parsephys( fea );
fea = parseprob( fea );

% Define manual source term with rank 3 to each vector compnent.
switch( opt.icase )
  case 1
    f1 = '0';
    f2 = '0';
    f3 = '(2*pi^2 + 1)*sin(pi*x)*sin(pi*y)';
end
fea.eqn.f.form{1} = 1;
fea.eqn.f.coef{1} = {cat(3,{f1},{f2},{f3})};


% Solve problem.
fea.sol.u = solvestat( fea, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot )
  postplot( fea, 'sliceexpr', 'Ec#1' )
  title( 'Curl solution #1')
end


% Error checking.
switch( opt.icase )
  case 1
    refsol = { '0', ...
               '0', ...
               'sin(pi*x)*sin(pi*y)', ...
               'pi*cos(pi*y)*sin(pi*x)', ...
               '-pi*cos(pi*x)*sin(pi*y)', ...
               '0' };
end
out.err = [ intsubd(['(',refsol{1},'-E#1)^2'],fea), ...
            intsubd(['(',refsol{2},'-E#2)^2'],fea), ...
            intsubd(['(',refsol{3},'-E#3)^2'],fea), ...
            intsubd(['(',refsol{4},'-Ec#1)^2'],fea), ...
            intsubd(['(',refsol{5},'-Ec#2)^2'],fea), ...
            intsubd(['(',refsol{6},'-Ec#3)^2'],fea) ];
out.pass = all( out.err < opt.tol );
if( nargout==0 )
  clear fea out
end
