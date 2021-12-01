function [ fea, out ] = ex_linearelasticity4( varargin )
%EX_LINEARELASTICITY4 Stress calculation of an I-beam attached to two brackets.
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY4( VARARGIN ) Example to calculate displacements and
%   stresses for an I-beam suppored by two brackets with circular holes.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {200e9}         Modulus of elasticity
%       nu          scalar {0.3}           Poissons ratio
%       force       scalar {1e5}           Load force
%       l           scalar {0.4}           Length of I-beam
%       ilev        scalar {2}             Grid regfinement level
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        200e9; ...
  'nu',       0.3; ...
  'force',    1e5; ...
  'l',        0.4; ...
  'ilev',     2; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      0.42; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
fea.sdim = { 'x' 'y' 'z' };   % Coordinate names.


% Grid generation.
fea.grid = get_grid( opt.ilev );


% Problem definition.
fea = addphys(fea,@linearelasticity);
fea.phys.el.eqn.coef{1,end} = { opt.nu };
fea.phys.el.eqn.coef{2,end} = { opt.E  };
fea.phys.el.sfun            = { opt.sfun opt.sfun opt.sfun };


% Boundary conditions.
dtol     = sqrt(eps);
fixbdr   = findbdr( fea, ['(sqrt(x.^2+z.^2)<=0.03+sqrt(eps))&(z>=-sqrt(eps))'] );
forcebdr = findbdr( fea, ['abs(y)>=0.2-sqrt(eps)'] );


% Fix boundaries (set zero Dirichlet BCs).
n_bdr  = max(fea.grid.b(3,:));        % Number of boundaries.
bctype = num2cell( zeros(3,n_bdr) );  % First set homogenous Neumann BCs everywhere.
[bctype{:,fixbdr}] = deal( 1 );       % Set Dirchlet BCs for right boundary.
fea.phys.el.bdr.coef{1,5} = bctype;

% Apply negative z-load to left boundary.
bccoef = num2cell( zeros(3,n_bdr) );
[bccoef{3,forcebdr}] = deal(-opt.force);
fea.phys.el.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );
fea.sol.u = solvestat( fea, 'fid', fid );


% Postprocessing.
if ( opt.iplot>0 )
  DSCALE = 5000;

  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2+w^2)', 'linestyle', 'none' )
  title( 'Total displacement' )
  view([30 20])

  subplot(1,2,2)
  dp = zeros(size(fea.grid.p));
  for i=1:3
    dp(i,:) = DSCALE*evalexpr( fea.dvar{i}, fea.grid.p, fea );
  end
  fea_disp.grid   = fea.grid;
  fea_disp.grid.p = fea_disp.grid.p + dp;
  plotgrid( fea_disp )
  title(['Displacement plot (at ',num2str(DSCALE),' times scale)'])
  view([30 20])

end


% Error check.
disp_max_ref = 6.204e-6;
xdisp = fea.sol.u(fea.eqn.dofm{1}(:));
ydisp = fea.sol.u(fea.eqn.dofm{2}(:)+fea.eqn.ndof(1));
zdisp = fea.sol.u(fea.eqn.dofm{3}(:)+sum(fea.eqn.ndof(1:2)));
disp  = sqrt(xdisp.^2+ydisp.^2+zdisp.^2);
disp_max = max(disp);

svm_max_ref = 4.410e6;
svm = evalexpr( fea.phys.el.eqn.vars{1,2}, fea.grid.p, fea );
svm_max = max(svm);

out.disp_max = disp_max;
out.svm_max  = svm_max;
out.err(1)   = abs(disp_max - disp_max_ref)/abs(disp_max_ref);
out.err(2)   = abs(svm_max - svm_max_ref)/abs(svm_max_ref);
out.pass     = all(out.err<opt.tol);


if ( nargout==0 )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ grid ] = get_grid( ilev )

n0 = 12;
nr = n0*2^(ilev-1);
r  = 0.03;   % Radius of bracket holes.
t  = 0.03;   % Thickness of brackets.
grid01 = ringgrid( 3*2^(ilev-1), 4*nr, r, r+t, [0;0] );
indc01 = selcells( grid01, 'y<=sqrt(eps)' );
grid01 = delcells( grid01, indc01 );

grid02 = holegrid( nr, 3*2^(ilev-1), (r+t)*[-1 1;-1 1], r, [0;0] );
indc02 = selcells( grid02, 'y>=-sqrt(eps)' );
grid02 = delcells( grid02, indc02 );
grid2d = gridmerge( grid01, findbdr(grid01,'y<=sqrt(eps)'), grid02, findbdr(grid02,'y>=-sqrt(eps)') );
t_br = 0.02;   % Width/depth of brackets.
d_br = 0.05;   % Separation distance between brackets.
grid1 = gridextrude( grid2d, 2^(ilev-1), t_br );
grid1 = gridrotate( grid1, pi/2, 1 );
grid2 = grid1;
grid1.p(2,:) = grid1.p(2,:) - d_br/2;
grid2.p(2,:) = grid2.p(2,:) + t_br + d_br/2;


% Create grids for the I-beam.
w_ib = 0.16;   % Beam width.
l_ib = 0.4;    % Beam length.
t_ib = 0.01;   % Beam thickness.
h_ib = 0.1;    % Beam height.

x_in    = linspace(-(r+t),r+t,n0+1);
x_coord = [ -w_ib/2 x_in w_ib/2];
y_coord = [ -0.2 -0.175 -0.15 -0.125 -0.1 -0.075 -(d_br/2+t_br) -d_br/2 0 d_br/2 d_br/2+t_br 0.075 0.1 0.125 0.15 0.175 0.2 ];
for i=2:ilev
  x_coord = sort([ x_coord [x_coord(1:end-1) + x_coord(2:end)]/2 ]);
  y_coord = sort([ y_coord [y_coord(1:end-1) + y_coord(2:end)]/2 ]);
end
grid3 = blockgrid( x_coord, y_coord, 2^(ilev-1), ...
                   [-w_ib/2 w_ib/2;-0.2 0.2;-(r+t)-t_ib -(r+t)] );
tm_ib = 2*(x_in(2)-x_in(1));
grid4 = blockgrid( 2*2^(ilev-1), y_coord, 5*2^(ilev-1), ...
                   [-tm_ib/2 tm_ib/2;-0.2 0.2;-(r+t)-t_ib-h_ib -(r+t)-t_ib] );
grid5 = grid3;
grid5.p(3,:) = grid5.p(3,:) - t_ib - h_ib;


% Merge grids.
tol = sqrt(eps)*1e3;
grid = gridmerge( grid1, findbdr(grid1,['z<=',num2str(-(r+t-tol))]), ...
                  grid3, findbdr(grid3,['z>=',num2str(-(r+t+tol))]) );
grid = gridmerge( grid2, findbdr(grid1,['z<=',num2str(-(r+t-tol))]), ...
                  grid,  findbdr(grid, ['(z>=',num2str(-(r+t+tol)),')&'...
                                        '(z<=',num2str(-(r+t-tol)),')']) );
grid = gridmerge( grid,  findbdr(grid,  ['z<=',num2str(-(r+t+t_ib-tol))]), ...
                  grid4, findbdr(grid4, ['z>=',num2str(-(r+t+t_ib+tol))]), 1 );
grid = gridmerge( grid,  findbdr(grid,  ['z<=',num2str(-(r+t+t_ib+h_ib-tol))]), ...
                  grid5, findbdr(grid5, ['z>=',num2str(-(r+t+t_ib+h_ib+tol))]), 2 );
