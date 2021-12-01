function [ fea, out ] = ex_periodic2( varargin )
%EX_PERIODIC2 2D Periodic Poisson equation example.
%
%   [ FEA, OUT ] = EX_PERIODIC2( VARARGIN ) 2D Periodic Poisson equation example.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       isol        scalar 0/{1}           Stationary/Time dependent solution
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'isol',     0;
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Unit square grid.
fea.sdim = {'x' 'y'};
fea.grid = rectgrid( 40 );


% Define Poisson equation with custom physics mode.
fea = addphys( fea, @customeqn, { 'u' } );
fea.phys.ce.eqn.seqn = '- ux_x - uy_y = f + u*eps';


% Boundary conditions.
DirBC = true;
NeuBC = ~DirBC;
NeuBC_coef = @periodic_bc;   % Assign periodic_bc solve hook function as boundary coefficient.
fea.phys.ce.bdr.coef = { 'bcnd_ce', '', '', {}, ...
                         { DirBC, NeuBC,      DirBC, NeuBC }, [], ...
                         { 0    , NeuBC_coef, 0,     0     } };

% Source term.
fea.expr = { 'f' 'x*sin(5*pi*y) + exp(-((x-0.5)^2 + (y-0.5)^2)/0.02)' };


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
if( opt.isol==0 )
  fea.sol.u = solvestat( fea, 'fid', opt.fid );
else
  fea.sol.u = solvetime( fea, 'tstop', -eps, 'fid', opt.fid );
end


% Postprocessing.
y  = linspace( 0, 1, 100 );
ul = evalexpr( 'u', [zeros(1,100);y], fea );
ur = evalexpr( 'u', [ones(1,100);y],  fea );
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot(fea, 'surfexpr', 'u', 'surfhexpr', 'u' )

  % Plot solution on left and right sides.
  hold on
  fea.grid.p(1,:) = fea.grid.p(1,:) - 1;
  postplot(fea, 'surfexpr', 'u', 'surfhexpr', 'u' )
  fea.grid.p(1,:) = fea.grid.p(1,:) + 2;
  postplot(fea, 'surfexpr', 'u', 'surfhexpr', 'u' )

  subplot(1,2,2)
  plot( y, ul, 'k-' )
  hold on
  plot( y, ur, 'r.' )
  grid on
  xlabel( 'y' )
  ylabel( 'u' )
  legend( 'u(x=0)', 'u(x=1)', 'Location', 'South' )
end


% Error checking.
err = norm( ur - ul );

out.err  = err;
out.pass = err<1e-6;
if ( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ A, f, prob ] = periodic_bc( A, f, prob, i_dvar, j_bdr, solve_step )

if( solve_step~=1 )   % Only process directly before linear solver.
  return
end
if( isstruct(A) )     % Sparsify matrix if necessary.
  A = sparse( A.rows, A.cols, A.vals, A.size(1), A.size(2) );
end


bdrm = prob.bdr.bdrm{i_dvar};
if( ~isempty(bdrm) )

  ndof = prob.eqn.ndof;
  if( ~iscell(ndof) )
    ndof = { ndof };
  end
  n_dof = sum( [ndof{1:(i_dvar-1)}] );

  % Find dofs and corresponding coordinates on boundary j.
  i_bdr    = mod(j_bdr+1,4) + 1;       % Hard wire i_bdr to j_bdr+2.
  ind      = find(bdrm(3,:)==i_bdr);
  [~,indh] = unique(bdrm(4,ind));      % Remove duplicate/shared points (optional).
  ind      = sort( ind(indh) );        % Unique and sorted index.
  indc     = bdrm(1,ind);
  inde     = bdrm(2,ind);
  dof_i_list = bdrm(4,ind) + n_dof;
  xii      = bdrm(6:end,ind);
  for k=1:numel(ind)
    y_i(k) = evalexpr0( 'y', xii(:,k)', [], indc(k), inde(k), prob );
  end

  % Build list of Dirchlet BC dofs and set right hand side/load vector to Dirchlet values.

  % Index to dofs with Dirichlet BCs on boundary j_bdr.
  ind      = find(bdrm(3,:)==j_bdr);   % Index to dofs on boundary j_bdr.
  [~,indh] = unique(bdrm(4,ind));      % Remove duplicate/shared points (optional).
  ind      = sort( ind(indh) );        % Unique and sorted index.
  for j=1:numel(ind)

    ix = ind(j);
    indc  = bdrm(1,ix);
    inde  = bdrm(2,ix);
    dof_j = bdrm(4,ix) + n_dof;
    xi = bdrm(6:end,ix);

    y = evalexpr0( 'y', xi', [], indc, inde, prob );
    [~,ix_i] = min( abs(y_i(:)-y) );
    dof_i = dof_i_list(ix_i);

    % Add i and j rows (to preserve both i and j equations).
    A(dof_i,:) = A(dof_i,:) + A(dof_j,:);
    f(dof_i)   = f(dof_i)   + f(dof_j);

    % Modify j rows so that A(j,i) - A(j,j) = 0.
    A(dof_j,:) = 0;   % Zero out periodic bc rows.
    A(dof_j,[dof_i dof_j]) =  [ 1 -1 ];
    f(dof_j) = 0;
  end
end
