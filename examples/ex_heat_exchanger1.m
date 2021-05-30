function [ fea, out ] = ex_heat_exchanger1( varargin )
%EX_FREE_CONVECTION Example of free and forced convection around a heated cylinder.
%
%   [ FEA, OUT ] = EX_HEAT_EXCHANGER1( VARARGIN ) Sets up and solves an example
%   of temperature transport in a fluid around a heated cylinder. The model
%   involves both free and forced convection by using the Boussinesq approximation.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.0005}        Max grid cell size
%       sf_u        string {sflag2}        Shape function for velocity
%       sf_p        string {sflag1}        Shape function for pressure
%       sf_T        string {sflag1}        Shape function for temperature
%       iplot       scalar 0/{1}           Plot solution and error (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.

cOptDef = { 'hmax',     0.0005;
            'sf_u',     'sflag1';
            'sf_p',     'sflag1';
            'sf_T',     'sflag1';
            'iplot',    1;
            'w',        0.0075;
            'h',        0.05;
            'yc',       0.02;
            'r',        0.003;
            'rho',      2.2e1;
            'miu',      2.8e-3;
            'cp',       3.1e3;
            'k',        0.55;
            'al',       0.26e-3;
            'g',        9.81;
            'vin',      40e-2;
            'Tc',       300;
            'Th',       330;
            'tol',      0.1;
            'fid',      1 };
[got,opt] = parseopt( cOptDef, varargin{:} );
fid       = opt.fid;


% Geometry definition.
fea.sdim         = { 'x' 'y' };
fea.geom.objects = { gobj_rectangle( 0, opt.w, 0, opt.h, 'R1' ) ...
                     gobj_circle( [0 opt.yc], opt.r, 'C1' ) };
fea = geom_apply_formula( fea, 'R1-C1' );


% Grid generation.
fea.grid = gridgen( fea, 'hmax', opt.hmax, 'fid', fid );


% Boundary conditions.
dtol = 1e-6;
ib_l = findbdr( fea, ['x<',num2str(dtol)] );         % Right boundary numbers.
ib_r = findbdr( fea, ['x>',num2str(opt.w-dtol)] );   % Left boundary number.
ib_b = findbdr( fea, ['y<',num2str(dtol)] );         % Bottom boundary number.
ib_t = findbdr( fea, ['y>',num2str(opt.h-dtol)] );   % Top boundary number.
ib_c = findbdr( fea, ['sqrt(x.^2+(y-',num2str(opt.yc),').^2)<',num2str(opt.r+dtol)] );   % Cylinder boundary numbers.


% Problem definition.
fea.expr = { 'rho'   opt.rho ;
             'miu'   opt.miu ;
             'cp'    opt.cp  ;
             'k'     opt.k   ;
             'alpha' opt.al  ;
             'g'     opt.g   ;
             'vin'   opt.vin ;
             'Tc'    opt.Tc  ;
             'Th'    opt.Th  };

fea = addphys(fea,@navierstokes);     % Add Navier-Stokes equations physics mode.
fea.phys.ns.eqn.coef{1,end} = { 'rho' };
fea.phys.ns.eqn.coef{2,end} = { 'miu' };
fea.phys.ns.eqn.coef{4,end} = { 'alpha*g*rho*(T-Tc)' };
fea.phys.ns.bdr.sel(ib_t)   = 4;      % Set pressure to zero on top boundary segment.
fea.phys.ns.bdr.sel(ib_b)   = 2;      % Set y-velocity to vin at bottom boundary.
fea.phys.ns.bdr.coef{2,end}{2,ib_b} = 'vin';
fea.phys.ns.sfun            = { opt.sf_u opt.sf_u opt.sf_p };

fea = addphys(fea,@heattransfer);     % Add heat transfer physics mode.
fea.phys.ht.sfun            = { opt.sf_T };
fea.phys.ht.eqn.coef{1,end} = { 'rho' };
fea.phys.ht.eqn.coef{2,end} = { 'cp' };
fea.phys.ht.eqn.coef{3,end} = { 'k' };
fea.phys.ht.eqn.coef{4,end} = { fea.phys.ns.dvar{1} };
fea.phys.ht.eqn.coef{5,end} = { fea.phys.ns.dvar{2} };
fea.phys.ht.bdr.sel         = 3*ones( 1, max(fea.grid.b(3,:)) );   % Default to symmetry BCs.
fea.phys.ht.bdr.sel(ib_t)   = 2;   % Set top boundary to convective flux/outflow.
fea.phys.ht.bdr.sel([ib_c ib_b])    = 1;   % Set temperature to bottom and cylinder boundaries.
[fea.phys.ht.bdr.coef{1,end}{ib_b}] = deal('Tc');
[fea.phys.ht.bdr.coef{1,end}{ib_c}] = deal('Th');


% Parse and solve problem.
fea = parsephys(fea);
[fea.bdr.d{1}{[ib_l ib_r]}] = deal(0);    % Manually apply slip boundary conditions.
[fea.bdr.d{2}{[ib_l ib_r]}] = deal([]);
[fea.bdr.n{2}{[ib_l ib_r]}] = deal(0);
fea = parseprob(fea);
fea.sol.u = solvestat( fea, 'fid', fid, 'maxnit', 50 );


% Postprocessing.
if( opt.iplot>0 )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'} )
  hold on
  fea.grid.p(1,:) = -fea.grid.p(1,:);   % Mirror solution.
  postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u' 'v'} )
  title('Velocity field')

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'T' )
  hold on
  fea.grid.p(1,:) = -fea.grid.p(1,:);   % Mirror solution.
  postplot( fea, 'surfexpr', 'T' )
  title('Temperature')
end

% Average temperature at outlet.
out.T_av_out = intbdr( 'T*v', fea, ib_t )/intbdr( 'v', fea, ib_t );
out.err = abs(intbdr('T-Tc',fea,3)/0.0075-1.7)/1.7;
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end
