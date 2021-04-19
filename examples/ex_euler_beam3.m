function [ fea, out ] = ex_euler_beam3( varargin )
%EX_EULER_BEAM3 1D Euler-Bernoulli beam vibration example.
%
%   [ FEA, OUT ] = EX_EULER_BEAM3( VARARGIN ) 1D Euler-Bernoulli beam model example for.
%   calculating vibration modes and eigenfrequencies. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       L           scalar {1}             Beam length
%       E           scalar {1}             Elastic modulus
%       I           expression {1}         Cross section moment of intertia
%       rho         expression {1}         Material density
%       A           expression {1}         Cross section area
%       nx          scalar {100}           Number of grid cells
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
            'rho',      1;
            'A',        6;
            'nx'        100;
            'iplot',    1;
            'tol',      1e-3;
            'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Grid generation.
fea.sdim = {'x'};
fea.grid = linegrid( opt.nx, 0, opt.L );


% Problem and equation definitions.
% Problem and equation definitions.
fea = addphys( fea, @eulerbeam );
fea.phys.eb.eqn.coef{1,end}  = { opt.rho };
fea.phys.eb.eqn.coef{2,end}  = { opt.A };
fea.phys.eb.eqn.coef{3,end}  = { opt.E };
fea.phys.eb.eqn.coef{4,end}  = { opt.I };
fea.phys.eb.bdr.coef{1,5}{1} = 1;


% Coefficients and equation/postprocessing expressions.
fea.expr = { 'L',    opt.L ;
             'E',    opt.E ;
             'I',    opt.I ;
             'rho',  opt.rho ;
             'A',    opt.A };


% Solve for eigenvalues and eigenvectors.
fea = parsephys( fea );
fea = parseprob( fea );

[fea.sol.u,fea.sol.l] = solveeig( fea, 'fid', opt.fid );
efq = sqrt(fea.sol.l)/(2*pi);


% Postprocessing.
if( opt.iplot )
  n_modes = 4;
  figure, hold on, grid on
  cols = {'b' 'r' 'y' 'g' 'c' 'm'};
  for i=1:n_modes
    postplot( fea, 'surfexpr', 'v', 'solnum', i, 'color', cols{mod(i-1,6)+1}, ...
              'axequal', 'off', 'grid', 'on', 'linewidth', 3 )
  end
  title( ['Vibration modes 1 - ',num2str(n_modes)] )
  xlabel( 'x' )
  ylabel( 'Displacement' )
end

% Error checking.
efq_ref = (pi*[0.59686,1.49418,2.50025,3.49999].').^2* ...
          sqrt(opt.E*opt.I/(opt.rho*opt.A*opt.L^4))/(2*pi);
err = abs(sort(efq(1:length(efq_ref)))-efq_ref)./efq_ref;


out.efq  = sort(efq);
out.err  = err;
out.pass = all(out.err<opt.tol);
if ( nargout==0 )
  clear fea out
end
