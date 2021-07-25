function [ fea, out ] = ex_resistive_heating1( varargin )
%EX_RESISTIVE_HEATING1 Resistive heating of a tungsten filament.
%
%   [ FEA, OUT ] = EX_RESISTIVE_HEATING1( VARARGIN ) Example to calculate temperature generated
%   by resistive heating of a tungsten spiral filament, such as for example in an incandescent
%   bulb. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       ibc         scalar 1{2}            Boundary type (1=insulation, 2=radiation)
%       isolve      scalar 1{2}            Solve type (1=stat.coupled, 2=split quasistatic)
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'ibc',      2;
            'isolve',   2;
            'sfun',     'sflag1';
            'iplot',    1;
            'tol',      1e-2;
            'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Set up geometry and grid.
r_wire   = 0.0005;
l_spiral = 0.005;
r_spiral = 0.004;

fea.sdim = { 'x' 'y' 'z' };
fea.grid = gridrotate( gridrevolve( circgrid( 2, 3, r_wire ), 80, l_spiral, 3, r_spiral ), -pi/2, 2 );
n_bdr = max(fea.grid.b(3,:));
ib2 = findbdr( fea, 'y<1e-3 & y>-1e-3' );
ib1 = setdiff(1:n_bdr,ib2);


% Problem constants and definitions.
V0    = 0.2;         % Voltage difference.
sigma = 1/52.8e-9;   % Electrical conductivity [1/Ohm/m].

rho   = 19.25e3;     % Density [kg/m3].
cp    = 133.9776;    % Specific heat capacity [J/kg/K].
k     = 173;         % Thermal conductivity [W/m/K].
if( opt.isolve==1 )
  Q   = [num2str(sigma),'*(Vx^2+Vy^2+Vz^2)'];   % Resistive heating source term.
else
  Q   = 'qfunction'; % Resistive heating source term function.
end

Const = 5.670367e-8;       % Constant for radiation BC (Stefan-Bolzmann constant).
T0    = 20+273.15;         % Initial temperature.
Tamb  = 80+273.15;         % Mean ambient temperature.

icub  = str2num(opt.sfun(end))+1; % Numerical quadrature order.
tmax  = 20;                 % Maximum time for temperature simulation.
dt    = tmax/20;            % Time step size.


% Set up and solve electrostatics problem.
feaV = fea;
feaV = addphys( feaV, @conductivemediadc );
feaV.phys.dc.sfun = { opt.sfun };

feaV.phys.dc.eqn.coef{2,end} = {sigma};    % Electrical conductivity.

feaV.phys.dc.bdr.sel(ib2) = 1;                 % Use Dirichlet / Electric potential BCs for end boundaries 57 and 58.
feaV.phys.dc.bdr.coef{1,end}{ib2(1)} = V0;     % Set electric potential on boundary 1 to V0 (boundary 57 is 0 by default).

if( opt.isolve==2 )
  feaV = parsephys( feaV );
  feaV = parseprob( feaV );
  u_V  = solvestat( feaV, 'fid', opt.fid );
end

% Set up and solve heat transfer problem.
fea = addphys( fea, @heattransfer );
fea.phys.ht.sfun = { opt.sfun };

fea.phys.ht.eqn.coef{1,end} = {rho};       % Density.
fea.phys.ht.eqn.coef{2,end} = {cp};        % Heat capacity.
fea.phys.ht.eqn.coef{3,end} = {k};         % Thermal conductivity.
fea.phys.ht.eqn.coef{7,end} = {Q};         % Heat source term.

if( opt.ibc==1 )
  fea.phys.ht.bdr.sel(ib1) = deal(3);     % Use insulation flux boundary condition for all boundaries except the ends.
else
  fea.phys.ht.bdr.sel(ib1) = deal(4);     % Use generalized heat flux boundary condition for all boundaries except the ends.
[fea.phys.ht.bdr.coef{4,end}{ib1}] = deal({0 0 0 Const Tamb});
end
fea.phys.ht.bdr.sel(ib2) = deal(1);     % Prescribe fixed ambient temperature at the end points.
fea.phys.ht.bdr.coef{1,end}{ib2(1)} = T0;
fea.phys.ht.bdr.coef{1,end}{ib2(2)} = T0;

if( opt.isolve==2 )
  % Add sigma and electric potential solution used by qfunction.
  fea.expr = { 's_sigma' '' '' sigma ;
               'u_V'     '' '' u_V   };

  fea  = parsephys( fea );
  fea  = parseprob( fea );
  l_write_qfunction()
  u_T  = solvetime( fea, 'icub', icub, 'init', {T0}, 'tstep', dt, 'tmax', tmax, 'fid', opt.fid );
  delete( fullfile(tempdir(),'qfunction.m') );
  clear qfunction;
end


% Merge solutions.
fea.phys.dc = feaV.phys.dc;
if( isfield(fea,'dvar') )
  fea = rmfield( fea, 'dvar' );
  fea = rmfield( fea, 'eqn'  );
end
fea = parsephys( fea );
fea = parseprob( fea );
if( opt.isolve==1 )
  fea.sol.u = solvestat( fea, 'icub', icub, 'fid', opt.fid, 'nlrlx', 1-0.5*(opt.ibc==2) );
  u_T = fea.sol.u(1:fea.eqn.ndof(1),:);
else
  fea.sol.u = [ u_T; repmat(u_V,1,size(u_T,2)) ];
end


% Postprocessing.
if( opt.iplot>0 )
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'V' )
  title( 'Electric potential, V')

  subplot(1,2,2)
  postplot( fea, 'surfexpr', 'T-273.15' )
  title( 'Temperature, T [C]')
end


% Error checking.
if( opt.isolve==1 )
  if( opt.ibc==1 )
    T_max_ref = 840;
  else
    T_max_ref = 688;
  end
else
  if( opt.ibc==1 )
    T_max_ref = 787;
  else
    T_max_ref = 680;
  end
end
out.err  = abs( max(u_T(:)) - T_max_ref )/T_max_ref;
out.pass = out.err < opt.tol;


if( nargout==0 )
  clear fea out
end



%-----------------------------
function [ Q ] = l_qfunction()
% Build up a source term expression from the precomputed solution
% vector stored in prob.expr{ u_V }. Expected to be called from evalexpr0.

Q = 'qfunction';

try
  % Extract variable from the calling function, evalexpr0.
  xi    = evalin( 'caller', 'xi' );
  ind_s = evalin( 'caller', 'ind_s' );
  ind_c = evalin( 'caller', 'ind_c' );
  prob  = evalin( 'caller', 'prob' );
  aJac  = evalin( 'caller', 'aJac' );

  n_sdim = size( prob.grid.p, 1 );
  n_vert = size( prob.grid.c, 1 );


  i_dvar = 1;   % Assume same dep. var. fem basis function for u_V as for i_dvar=1.
  [~,~,~,sfun] = evalsfun( prob.sfun{i_dvar}, 0, n_sdim, n_vert );   % Get shape function root string.
  if( isempty(aJac) )   % Calculate Jacobian if required.
    store_aJTmp = strcmpi(sfun(end-1:end),'H3') && ( (n_sdim==2 && n_vert==4) || (n_sdim==3 && n_vert==8) );
    [~,aJac] = tfjac( 1, prob.grid.p, prob.grid.c(:,ind_c), [], xi, aJac, store_aJTmp );
  end

  u_V = prob.expr{ find(strcmp(prob.expr,'u_V')), end };   % V solution vector.

  % Evaluate derivatives of V (1=function value, 2=x-derivative, 3=y-derivative, 4=z-derivative).
  Vx  = evaldvar( sfun, 2, n_sdim, n_vert, xi, aJac, prob.eqn.dofm{i_dvar}(:,ind_c), u_V );
  Vy  = evaldvar( sfun, 3, n_sdim, n_vert, xi, aJac, prob.eqn.dofm{i_dvar}(:,ind_c), u_V );
  Vz  = evaldvar( sfun, 4, n_sdim, n_vert, xi, aJac, prob.eqn.dofm{i_dvar}(:,ind_c), u_V );

  % Expression for sigma.
  s_sigma = prob.expr{ find(strcmp(prob.expr,'s_sigma')), end };

  % Calculate final expression and return values.
  Q = evalexpr0( s_sigma, xi, ind_s, ind_c, [], prob, aJac ).*( Vx.^2 + Vy.^2 + Vz.^2 );
catch,end

%-----------------------------
function l_write_qfunction()

ds    = dbstack;
file  = which(ds(1).file);
s     = fileread( file );
poss  = strfind( s, 'function [ Q ] = l_qfunction' );
pose  = strfind( s, 'function l_write_qfunction' );
s     = strrep( s(poss:pose-1), 'l_qfunction', 'qfunction' );
fid   = fopen( fullfile(tempdir(),'qfunction.m'), 'w' );
fprintf( fid, '%s', s );
fclose( fid );
addpath( tempdir() );
