function [ fea, out ] = ex_heattransfer10( varargin )
%EX_HEATTRANSFER10 Conjugate heat transfer test example with multiple domains
%
%   [ FEA, OUT ] = EX_HEATTRANSFER10( VARARGIN ) Sets up and solves a
%   multiple domain conjugate heat transfer example.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.0005}        Max grid cell size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       sf_T        string {sflag1}        Shape function for temperature
%       solver      string openfoam/{}     Use OpenFOAM or default solver
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2026 Precise Simulation, Ltd.

cOptDef = { 'hmax',     1/10;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'sf_T',     'sflag1';
            'turb',     [];
            'solver',   '';
            'iplot',    1;
            'tol',      0.1;
            'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
fea.sdim         = { 'x' 'y' };
fea.geom.objects = { gobj_rectangle(0, 1, 0, 0.5), ...
                     gobj_rectangle(0, 1, 0.5, 1, 'R2') };


% Grid generation.
n = round(1/opt.hmax);
fea.grid = gridmerge( rectgrid(2*n, n, [0, 1; 0, 0.5]), 3, rectgrid(2*n, n, [0, 1; 0.5, 1]), 1);
fea.grid.s( selcells(fea,'y >= 0.5') ) = 2;
fea.grid = gridbdrx(fea.grid);


% Problem definition.
fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
fea.phys.ns.eqn.coef{1,end} = {1, 1000};
fea.phys.ns.eqn.coef{2,end} = {1, 959e-6};
fea.phys.ns.prop.active(:,1) = 0;
fea.phys.ns.bdr.sel(6) = 2;
fea.phys.ns.bdr.coef{2,end}{1,6} = 0.00001;
fea.phys.ns.bdr.sel(4) = 4;

fea = addphys(fea,@heattransfer);     % Add heat transfer physics mode.
fea.phys.ht.eqn.coef{1,end} = {8000, 1000};
fea.phys.ht.eqn.coef{2,end} = {450, 4181};
fea.phys.ht.eqn.coef{3,end} = {80, 4181*959e-6/6.62};
fea.phys.ht.eqn.coef{4,end} = {0, 'u'};
fea.phys.ht.eqn.coef{5,end} = {0, 'v'};
fea.phys.ht.bdr.sel([1,6]) = 1;
fea.phys.ht.bdr.coef{1,end}{1} = 400;
fea.phys.ht.bdr.coef{1,end}{6} = 300;
fea.phys.ht.bdr.sel([2,3,5]) = 3;
fea.phys.ht.eqn.coef{end}{1} = 300;
fea.phys.ht.eqn.coef{end}{2} = 300;


% Parse and solve problem.
fea = parsephys(fea);
fea = parseprob(fea);

if( strcmp(opt.solver,'openfoam') )
  fea.sol.u = openfoam( fea, 'fid', fid, 'ddtScheme', 'CrankNicolson', 'endTime', 1e5, 'maxDeltaT', 1e5/100, 'nproc', 1, 'turb', opt.turb );
else
  fea.sol.u = solvestat( fea, 'fid', fid, 'maxnit', 50 );
end


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'} )
  title('Velocity field')

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'T' )
  title('Temperature')
end

% Average temperature at outlet.
out.ref = [1.4111e-5, 309.7998, 394.3849];
out.val = [evalexpr('u',[1;0.75],fea), ...
           evalexpr('T',[1;0.75],fea), ...
           evalexpr('T',[1;0.5],fea)];
out.err = abs(out.val - out.ref)./out.ref;
out.pass = all(out.err < opt.tol);

if( nargout==0 )
  clear fea out
end
