function [ fea, out ] = ex_euler_beam2( varargin )
%EX_EULER_BEAM2 1D Euler-Bernoulli beam model example.
%
%   [ FEA, OUT ] = EX_EULER_BEAM2( VARARGIN ) 1D Euler-Bernoulli beam model example
%   with a distributed load. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       L           scalar {2}             Beam length
%       E           scalar {3}             Elastic modulus
%       I           expression {4}         Cross section moment of intertia
%       q           expression {-5}        Beam force
%       nx          scalar {6}             Number of grid cells
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'L',        2;
            'E',        3;
            'I',        4;
            'q',       -5;
            'nx'        6;
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Grid generation.
fea.sdim = {'x'};
fea.grid = linegrid( opt.nx, 0, opt.L );


% Problem and equation definitions.
fea = addphys( fea, @eulerbeam );
fea.phys.eb.eqn.coef{3,end}  = { opt.E };
fea.phys.eb.eqn.coef{4,end}  = { opt.I };
fea.phys.eb.eqn.coef{5,end}  = { opt.q };
fea.phys.eb.bdr.coef{1,5}{1} = 1;
fea = parsephys( fea );


% Coefficients and equation/postprocessing expressions.
fea.expr = { 'L',  opt.L ;
             'M',  fea.phys.eb.eqn.vars{3,2} };


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 3, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot )
  figure
  subplot(3,1,1), hold on
  postplot( fea, 'surfexpr', 'E_eb*I_eb*v/(q_eb*L^4)', 'linewidth', 2 )
  postplot( fea, 'surfexpr', 'x^2*(6*L^2-4*L*x+x^2)/24/L^4', 'color', 'r', 'linestyle', ':' )
  title( 'E_eb*I_eb*v(x)/(q_eb*L^4)' )
  axis normal, grid on

  subplot(3,1,2), hold on
  postplot( fea, 'surfexpr', 'E_eb*I_eb*vx/(q_eb*L^3)', 'linewidth', 2 )
  postplot( fea, 'surfexpr', 'x*(3*L^2-3*L*x+x^2)/6/L^3', 'color', 'r', 'linestyle', ':' )
  title( 'E_eb*I_eb*theta(x)/(q_eb*L^3)' )
  axis normal, grid on

  subplot(3,1,3), hold on
  postplot( fea, 'surfexpr', 'M/(q_eb*L^2)', 'linewidth', 2 )
  postplot( fea, 'surfexpr', '-1/2*(L-x)^2/L^2', 'color', 'r', 'linestyle', ':' )
  axis normal, grid on
  title( 'M(x)/(q_eb*L^2)' )
end


% Error checking.
err_v  = evalexpr( 'abs(v-q_eb*x^2*(6*L^2-4*L*x+x^2)/(24*E_eb*I_eb))', ...
                   linspace(0,opt.L,3*opt.nx), fea );
err_th = evalexpr( 'abs(vx-q_eb*x*(3*L^2-3*L*x+x^2)/(6*E_eb*I_eb))', ...
                   linspace(0,opt.L,3*opt.nx), fea );
err_M  = evalexpr( 'abs(vxx-q_eb/2*(L-x)^2/(E_eb*I_eb))', ...
                   linspace(0,opt.L,3*opt.nx), fea );
err = norm([err_v;err_th;err_M]);


out.err  = err;
out.pass = out.err<opt.tol;
if ( nargout==0 )
  clear fea out
end
