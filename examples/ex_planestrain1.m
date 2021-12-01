function [ fea, out ] = ex_planestrain1( varargin )
%EX_PLANESTRAIN1 Plane strain analysis of a pressure vessel.
%
%   [ FEA, OUT ] = EX_PLANESTRAIN1( VARARGIN ) Benchmark example for plane strain
%   approximation of a pressure vessel (annular cross section with symmetry).
%
%   Reference.  B. J. Mac Donald, Practical Stress Analysis with Finite Elements (2nd Ed),
%   case study E on page 327, 2007.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       E           scalar {207e9}         Modulus of elasticity
%       nu          scalar {0.27}          Poissons ratio
%       sfun        string {sflag1}        Shape function for displacements
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'E',        207e9; ...
  'nu',       0.27; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'igrid',    1; ...
  'tol',      0.1; ...
  'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry and grid.
fea.sdim = { 'x' 'y' };   % Coordinate names.
fea.grid = ringgrid( 12, 216, 100e-3, 120e-3 );
fea.grid = delcells( fea.grid, selcells( fea.grid, '(x<=eps) | (y<=eps)') );
if( opt.igrid~=1 )
  fea.grid = quad2tri( fea.grid );
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Problem definition.
fea = addphys( fea, @planestrain );
fea.phys.psn.eqn.coef{1,end} = { opt.nu };
fea.phys.psn.eqn.coef{2,end} = { opt.E  };
fea.phys.psn.sfun            = { opt.sfun opt.sfun };


% Boundary conditions.
bctype = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bctype{1,4} = 1;
bctype{2,3} = 1;
fea.phys.psn.bdr.coef{1,5} = bctype;

bccoef = mat2cell( zeros(2,n_bdr), [1 1], ones(1,n_bdr) );
bccoef{1,1} = '-nx*1e4';
bccoef{2,1} = '-ny*1e4';
fea.phys.psn.bdr.coef{1,end} = bccoef;


% Parse and solve problem.
fea       = parsephys( fea );
fea       = parseprob( fea );                          % Check and parse problem struct.
fea.sol.u = solvestat( fea, 'fid', fid, 'icub', 1+str2num(strrep(opt.sfun,'sflag','')) );   % Call to stationary solver.


% Postprocessing.
s_disp = fea.phys.psn.eqn.vars{2,end};
if( opt.iplot>0 )
  figure
  postplot( fea, 'surfexpr', s_disp )
  title( 'Total displacement' )
end


% Error checking.
s_sx     = fea.phys.psn.eqn.vars{5,end};
s_sy     = fea.phys.psn.eqn.vars{6,end};
s_sxy    = fea.phys.psn.eqn.vars{8,end};
s_sp1    = fea.phys.psn.eqn.vars{9,end};
s_sp2    = fea.phys.psn.eqn.vars{10,end};
s_sp3    = fea.phys.psn.eqn.vars{11,end};
v_disp   = evalexpr( s_disp, [100e-3 120e-3-2*sqrt(eps);0 0]+sqrt(eps), fea )';
v_dref   = [2.64e-8 2.41e-8];
[v_sx(1),v_sx(2)] = minmaxsubd( s_sx, fea );
v_sxref  = [-10000 55454];
[v_sy(1),v_sy(2)] = minmaxsubd( s_sy, fea );
v_syref  = [-10000 55454];
[v_sxy(1),v_sxy(2)] = minmaxsubd( s_sxy, fea );
v_sxyref = [-32730 0];
[v_sp1(1),v_sp1(2)] = minmaxsubd( s_sp1, fea );
v_sp1ref = [4.5e4 55454];
[v_sp2(1),v_sp2(2)] = minmaxsubd( s_sp2, fea );
v_sp2ref = [1.227e4 1.227e4];
[v_sp3(1),v_sp3(2)] = minmaxsubd( s_sp3, fea );
v_sp3ref = [-1e4 0];
out.err  = [ abs([v_dref-v_disp])./v_dref ;
             abs([v_sxref-v_sx])./v_sxref ;
             abs([v_syref-v_sy])./v_syref ;
             abs([v_sxyref(1)-v_sxy(1)])./v_sxyref(1) 0 ;
             abs([v_sp1ref(2)-v_sp1(2)])./v_sp1ref(2) 0 ;
             abs([v_sp2ref-v_sp2])./v_sp2ref      ;
             abs([v_sp3ref(1)-v_sp3(1)])./v_sp3ref(1) 0 ];
out.pass = all( out.err(:) <= opt.tol );


if( nargout==0 )
  clear fea out
end
