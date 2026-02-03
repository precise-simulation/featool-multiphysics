function [ fea, out ] = ex_heat_exchanger2( varargin )
%EX_HEAT_EXCHANGER2 Multi-Simulation heat exchanger.
%

%   [ FEA, OUT ] = EX_HEAT_EXCHANGER2( VARARGIN ) This model
%   illustrates the multi-solver simulation process for a coupled tube
%   and fin heat exchanger.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.2}           Max grid cell size
%       solver      string {}              Solver selection default or openfoam
%       coupled     logical 1{0}           Use coupled or segregated solution process
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2026 Precise Simulation, Ltd.

cOptDef = { 'hmax',   0.2;
            'solver', '';
            'coupled', 0;
            'tolres',  1e-3;
            'iplot',   1;
            'rho',     1;
            'miu',     0.00526;
            'k',       3.76e-3;
            'uin',     1;
            'tol',     0.15;
            'fid',     1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
fea.sdim = { 'x', 'y', 'z' };
fea.geom.objects = { gobj_block( 0, 20, 0, 5, 0, 1, 'B1' ), ...
                     gobj_cylinder( [10, 0, 0], 2.5, 1, [0, 0, 1], 'C1' ) };
fea = geom_apply_formula( fea, 'B1-C1' );
fea = geom_add_gobj( fea, gobj_block( 5, 15, 0, 5, 0, 0.0625, 'B2' ) );
fea.geom = copy_geometry_object( 'B2', fea.geom, {[0, 0, 0.9375]} );
fea = geom_apply_formula( fea, 'CS1-B2-B3' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );


% Flow problem definition.
fea = addphys( fea, @navierstokes );
fea.phys.ns.eqn.coef{1,end} = { opt.rho };
fea.phys.ns.eqn.coef{2,end} = { opt.miu };
fea.phys.ns.bdr.sel(5) = 2;  % Inlet.
fea.phys.ns.bdr.coef{2,end}{1,5} = opt.uin;
fea.phys.ns.bdr.sel(11) = 4;  % Outlet.
fea.phys.ns.bdr.sel([1:2, 4, 6, 12, 15, 17]) = 5;  % Symmetry/slip.
fea.phys.ns.prop.artstab.id = 1;


% Parse and solve flow problem.
if( ~opt.coupled )
  fea = parsephys( fea );
  fea = parseprob( fea );
  if( strcmp(opt.solver,'openfoam') )
    logfid = fid; if( ~got.fid ), fid = []; end
    fea.sol.u = openfoam( fea, 'tolres', opt.tolres, 'fid', fid, 'logfid', logfid );
    fid = logfid;
  else
    fea.sol.u = solvestat( fea, 'fid', fid, 'nlrlx', 0.9, 'maxnit', 50 );
  end

  fea.phys.ns.prop.active = zeros(4,1);
end


% Heat transfer problem definition.
fea = addphys( fea, @heattransfer, { 'T' } );
fea.phys.ht.eqn.coef{1,end} = { opt.rho };
fea.phys.ht.eqn.coef{3,end} = { opt.k };
fea.phys.ht.eqn.coef{4,end} = { 'u' };
fea.phys.ht.eqn.coef{5,end} = { 'v' };
fea.phys.ht.eqn.coef{6,end} = { 'w' };
fea.phys.ht.prop.artstab.sd = 1;
fea.phys.ht.prop.artstab.sd_coef = 1.0;
fea.phys.ht.bdr.sel(5) = 1;  % Tc
fea.phys.ht.bdr.coef{1,end}{1,5} = 0;
fea.phys.ht.bdr.sel([3, 7:10, 13:14, 16]) = 1;  % Th
[fea.phys.ht.bdr.coef{1,end}{1,[3, 7:10, 13:14, 16]}] = deal(1);
fea.phys.ht.bdr.sel(11) = 2;  % Outflow.
fea.phys.ht.bdr.sel([1:2, 4, 6, 12, 15, 17]) = 3;  % Insulation.


% Parse and solve problem.
fea = parsephys( fea );
fea = parseprob( fea );
if( opt.coupled && strcmp(opt.solver,'openfoam') )
  logfid = fid; if( ~got.fid ), fid = []; end
  fea.sol.u = openfoam( fea, 'tolres', opt.tolres, 'fid', fid, 'logfid', logfid );
  fid = logfid;
else
  solvargs = {};
  if( ~opt.coupled )
    solvargs = {'init', 'restart', 'solcomp', [0, 0, 0, 0, 1]'};
  end
  fea.sol.u = solvestat( fea, 'fid', fid, solvargs{:} );
end


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'sliceexpr', 'T', 'slicex', 10, 'slicey', 2.5, 'slicez', 0.5, 'title', 'slice: Temperature, T' );
end

out.Tref = 3.8;
out.Tcalc = intbdr('T',fea,11);
out.err = abs(out.Tref - out.Tcalc) / out.Tref;
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
