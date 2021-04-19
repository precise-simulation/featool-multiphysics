function [ fea, out ] = ex_euler_beam1( varargin )
%EX_EULER_BEAM1 1D Euler-Bernoulli beam model example.
%
%   [ FEA, OUT ] = EX_EULER_BEAM1( VARARGIN ) 1D Euler-Bernoulli beam model example,
%   cantilever beam with point load. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       L           scalar {2}             Beam length
%       E           scalar {3}             Elastic modulus
%       I           expression {4}         Cross section moment of intertia
%       q           expression {-5}        Point load
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
fea.phys.eb.bdr.coef{1,5}{1} = 1;
fea.phys.eb.bdr.coef{1,7}{2}{1} = opt.q;
fea = parsephys( fea );


% Coefficients and equation/postprocessing expressions.
fea.expr = { 'L',  opt.L ;
             'M',  fea.phys.eb.eqn.vars{3,2} ;
             'P',  opt.q ;
             'v_ref', 'P*x^2*(3*L-x)/(6*E_eb*I_eb)' };


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvestat( fea, 'icub', 3, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot )
  figure
  subplot(2,1,1), hold on
  postplot( fea, 'surfexpr', 'v', 'linewidth', 2 )
  postplot( fea, 'surfexpr', 'v_ref', 'color', 'r', 'linestyle', ':' )
  title( 'v(x)' )
  axis normal, grid on

  subplot(2,1,2), hold on
  postplot( fea, 'surfexpr', 'M', 'linewidth', 2 )
  postplot( fea, 'surfexpr', '-P*(L-x)', 'color', 'r', 'linestyle', ':' )
  axis normal, grid on
  title( 'M(x)' )
end


% Error checking.
err_v  = evalexpr( 'abs(v-v_ref)', ...
                   linspace(0,opt.L,3*opt.nx), fea );
err = norm(err_v);


out.err  = err;
out.pass = out.err<opt.tol;
if ( nargout==0 )
  clear fea out
end
