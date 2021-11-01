function [ fea, out ] = ex_periodic1( varargin )
%EX_PERIODIC1 1D Example of pulse in a periodic line.
%
%   [ FEA, OUT ] = EX_PERIODIC1( VARARGIN ) Moving 1D pulse in a periodic domain.
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       init        string '(x>=0.4)*(x<=0.6)'   Initial shape of pulse
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'init',     '(x>=0.4)*(x<=0.6)';
            'iplot',    1;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;

fea.sdim = {'x'};
fea.grid = linegrid( 100 );

fea = addphys( fea, @convectiondiffusion );
fea.phys.cd.eqn.coef{3,end} = { 1 };
fea.phys.cd.eqn.coef{2,end} = { 0.002 };
fea = parsephys( fea );

fea.bdr.d{1} = { [], [] };
fea.bdr.n{1}{2} = @periodic_bc;   % Assign periodic_bc solve hook function as boundary coefficient.


% Parse and solve problem.
fea = parseprob( fea );
fea.sol.u = solvetime( fea, 'init', opt.init, 'tstep', 0.005, 'tmax', 1, 'fid', opt.fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'c', 'solnum', 1, 'color', 'b' )
  postplot( fea, 'surfexpr', 'c', 'solnum', floor(size(fea.sol.u,2)/4), 'color', 'c' )
  postplot( fea, 'surfexpr', 'c', 'solnum', floor(size(fea.sol.u,2)/2), 'color', 'm' )
  postplot( fea, 'surfexpr', 'c', 'solnum', floor(size(fea.sol.u,2)*3/4), 'color', 'g' )
  postplot( fea, 'surfexpr', 'c', 'solnum', size(fea.sol.u,2) )
  axis( [0 1 0 1] )
  grid on
end


% Error checking.
if( strcmp(opt.init,'(x>=0.4)*(x<=0.6)') )
  i_ref  = 0.2;
  i_calc = intsubd( 'c*(x>=0.4)*(x<=0.6)', fea );
  err = abs(i_calc-i_ref)/i_ref;

  if( ~isempty(fid) )
    fprintf(fid,'\nError: %f\n',err)
    fprintf(fid,'\n\n')
  end
end


out.err  = err;
out.pass = err<0.3;
if ( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ A, f, prob ] = periodic_bc( A, f, prob, i_dvar, j_bdr, solve_step )

if( j_bdr~=2 )        % Hard wired to only process boundary 2.
  return
end
if( solve_step~=1 )   % Only process directly before linear solver.
  return
end
if( isstruct(A) )     % Sparsify matrix if necessary.
  A = sparse( A.rows, A.cols, A.vals, A.size(1), A.size(2) );
end

% Degrees of freedom to couple.
dof_i = 1;
dof_j = prob.eqn.ndof;

% Add i and j rows (to preserve both i and j equations).
A(dof_i,:) = A(dof_i,:) + A(dof_j,:);
f(dof_i)   = f(dof_i)   + f(dof_j);

% Modify j rows so that A(j,i) - A(j,j) = 0.
A(dof_j,:) = 0;   % Zero out periodic bc rows.
A(dof_j,[dof_i dof_j]) =  [ 1 -1 ];
f(dof_j) = 0;
