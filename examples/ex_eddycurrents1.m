function [ fea, out ] = ex_eddycurrents1( varargin )
%EX_EDDYCURRENTS1 2D Eddy currents test example.
%
%   [ FEA, OUT ] = EX_EDDYCURRENTS1( VARARGIN ) 2D Eddy currents test
%   example for vector elements (Nedelec). Accepts the following
%   property/value pairs.
%
%   Reference:
%
%   [1] I. Anjam, J. Valdman, Fast MATLAB assembly of FEM matrices in
%   2D and 3D: Edge elements, Applied Mathematics and Computation, vol
%   267, 2015, pp 252-263, DOI:10.1016/j.amc.2015.03.105.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       icase       integer 2/{1}          Predefined test case
%       hmax        scalar {0.01}          Grid cell size
%       sfun        string {sf_simp_N1}    Vector shape function (Nedelec)
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'icase',    1;
            'hmax',     0.01;
            'sfun',     'sf_simp_N1';
            'iplot',    1;
            'tol',      0.01;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});


% Geometry and grid generation.
fea.sdim = {'x','y'};
fea.geom.objects = { gobj_rectangle() };
fea.grid = quad2tri( rectgrid(ceil(1/opt.hmax)) );


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
    f1 = '2+y*(1-y)';
    f2 = '2+x*(1-x)';
  case 2
    f1 = '(x>y)*(sin(2*pi*x) + 2*pi*cos(2*pi*x)*(x - y)) - (x>y)*(cos(y*(x - y)^2*(x - 1)^2)*(2*y*(x - 1)^2 - (2*x - 2*y)*(x - 1)^2 - (2*x - 2)*(x - y)^2 + y*(2*x - 2*y)*(2*x - 2)) + sin(y*(x - y)^2*(x - 1)^2)*((x - y)^2*(x - 1)^2 - y*(2*x - 2*y)*(x - 1)^2)*(y*(2*x - 2*y)*(x - 1)^2 + y*(2*x - 2)*(x - y)^2))';
    f2 = '(x>y)*(sin(y*(x - y)^2*(x - 1)^2) - sin(2*pi*x)) + (x>y)*sin(y*(x - y)^2*(x - 1)^2)*(y*(2*x - 2*y)*(x - 1)^2 + y*(2*x - 2)*(x - y)^2)^2 - (x>y)*cos(y*(x - y)^2*(x - 1)^2)*(2*y*(x - 1)^2 + 2*y*(x - y)^2 + 2*y*(2*x - 2*y)*(2*x - 2)) + 4*(x>y)*pi^2*sin(2*pi*x) - 4*(x>y)*pi^2*sin(2*pi*x)';
end
fea.eqn.f.form{1} = 1;
fea.eqn.f.coef{1} = {cat(3,{f1},{f2})};


% Solve problem.
fea.sol.u = solvestat( fea, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot )
  subplot(2,2,1)
  postplot( fea, 'surfexpr', 'E#1', 'surfhexpr', 'E#1' )
  title( 'Solution component #1')
  subplot(2,2,2)
  postplot( fea, 'surfexpr', 'E#2', 'surfhexpr', 'E#2' )
  title( 'Solution component #2')
  subplot(2,2,3)
  postplot( fea, 'surfexpr', 'Ec', 'surfhexpr', 'Ec' )
  title( 'Curl solution')
end


% Error checking.
switch( opt.icase )
  case 1
    refsol = { 'y*(1-y)', ...
               'x*(1-x)', ...
               '2*(y-x)' };
  case 2
    refsol = { '(x>y)*(sin(2*pi*x) + 2*pi*cos(2*pi*x)*(x - y))', ...
               '(x>y)*(sin(y*(x-y)^2*(x-1)^2) - sin(2*pi*x))', ...
               '(x>y)*(2*y*(x-y)*(x-1)*(2*x-y-1)*cos(y*(x-y)^2*(x-1)^2))' };
end
out.err = [ intsubd(['(',refsol{1},'-E#1)^2'],fea), ...
            intsubd(['(',refsol{2},'-E#2)^2'],fea), ...
            intsubd(['(',refsol{3},'-Ec)^2'],fea) ];
out.pass = all( out.err < opt.tol );
if( nargout==0 )
  clear fea out
end
