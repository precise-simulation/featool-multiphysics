function [ fea, out ] = ex_potential_flow1( varargin )
%EX_POTENTIAL_FLOW1 Potential flow around a NACA wing profile.
%
%   [ FEA, OUT ] = EX_POTENTIAL_FLOW1( VARARGIN ) Example to calculate
%   potential flow around a NACA airfoil with the stream-function
%   formulation. Using fminbnd to establish the stream-function value
%   on the body so as to enforce the Kutta condition and minimize the
%   velocity at the trailing edge.
%
%   Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       series      string {0012}          NACA 4-series wing profile descriptor
%       alfa        scalar {6}             Angle of attack
%       sfun        string {sflag1}        Shape function for stream-function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct
%
%   See also FMINBND

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { 'series',   '0012';
            'alfa',      6;
            'sfun',     'sflag1';
            'iplot',     1;
            'tol',       0.05;
            'fid',       1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;

series = opt.series;
alfa   = opt.alfa;
sfun   = opt.sfun;


func = @(psi) velocity_at_trailing_edge( psi, series, alfa, sfun );
args = {};
if( ~isempty(opt.fid) )
  args = { optimset('Display','iter','TolX',1e-3) };
end
[psi_min,U_min] = fminbnd( func, -10, 10, args{:} );
% [psi_min,U_min] = fminsearch( func, 0, args{:} );

fea = l_setup_fea( psi_min, series, alfa, sfun );
fea.sol.u = solvestat( fea, 'fid', [] );
err = l_postprocessing( opt.iplot, fea, series, alfa );
out.err  = err;
out.pass = err < opt.tol;

if( ~nargout )
  clear fea out
end


%------------------------------------------------------------------------------%
function [ U_T ] = velocity_at_trailing_edge( psi_body, series, alfa, sfun )
% Compute velocity at the trailing edge with a constant
% stream-function value psi_body applied to the wing.

fea = l_setup_fea( psi_body, series, alfa, sfun );

fea.sol.u = solvestat( fea, 'fid', [] );

p_T = fea.geom.objects{1}.children{2}.v(end,:).';
U_T = sqrt(evalexpr( 'psix^2+psiy^2', p_T, fea ));

%------------------------------------------------------------------------------%
function [ fea_out ] = l_setup_fea( psi_body, series, alfa, sfun )
% Define a fea struct for potential flow around a wing profile with
% the stream-function formulation

persistent fea

if( nargin<4 )
  sfun = 'sflag1';
end
if( nargin<3 )
  alfa = 0;
end
if( nargin<2 )
  series = '0012';
end
if( nargin<1 )
  psi_body = 0;
end


if( isempty(fea) )

  % Constants and expressions.
  fea.expr = { 'u',   'psiy';
               'v',  '-psix';
               'U',   'sqrt(u^2+v^2)';
               'uinf', 1;
               'vinf', 0;
               'cp',  '1-U^2/(uinf^2+vinf^2)' };

  % Geometry.
  fea.sdim = { 'x', 'y' };
  gobj1 = gobj_naca( series, [0 0], 1, alfa, 100, 'N1' );
  gobj2 = gobj_rectangle( -1, 2, -1.5, 1.5, 'R1' );
  fea.geom.objects = {gobj1 gobj2};
  fea.geom = geom_apply_formula( fea.geom, 'R1-N1' );

  % Grid.
  HMAX  = 0.3;    % Largest grid size.
  HMAXB = 0.01;   % Grid size around wing profile.
  fea.grid = gridgen( fea, 'gridgen', 'gridgen2d', ...
                      'hmax', HMAX, 'hmaxb', [ HMAX*ones(1,4), HMAXB*ones(1,2) ], ...
                      'nsm', 25, 'q', 0.65, 'blayer', [ 0.1, 5 ], 'fid', [] );

  % Physics mode and equation settings.
  fea = addphys( fea, @poisson, {'psi'} );
  fea.phys.poi.sfun = { sfun };
  fea.phys.poi.eqn.coef{3,end} = {0};

  % Boundary conditions.
  i_bdr_tb   = [3,1];
  i_bdr_lr   = [4,2];
  i_bdr_body = [5,6];
  fea.phys.poi.bdr.sel([i_bdr_tb,i_bdr_body]) = 1;   % Dirichlet BCs for top/bottom/body boundaries.
  fea.phys.poi.bdr.sel([i_bdr_lr]) = 2;              % Neumann BCs for left/right boundaries.

  dbc = cell(1,6);
  [dbc{i_bdr_tb}]   = deal('uinf*y');
  [dbc{i_bdr_body}] = deal(psi_body);
  fea.phys.poi.bdr.coef{1,end} = dbc;

else

  % Update stream function value on wing profile body.
  fea = rmfield(fea,'bdr');

  dbc = cell(1,6);
  i_bdr_tb   = [3,1];
  [dbc{i_bdr_tb}]   = deal('uinf*y');
  i_bdr_body = [5,6];
  [dbc{i_bdr_body}] = deal(psi_body);
  fea.phys.poi.bdr.coef{1,end} = dbc;

end

fea = parsephys( fea );
fea = parseprob( fea );
fea_out = fea;

%------------------------------------------------------------------------------%
function [ err ] = l_postprocessing( iplot, fea, series, alfa )
% Postprocessing and visualization.

if( iplot )
  figure
  subplot(1,2,1)
  postplot( fea, 'surfexpr', 'U', 'isoexpr', 'psi', 'isolev', 50, ...
            'colorbar', false, 'axis', false )
  title( 'Streamfunction \psi, and velocity field \it{U}' )
end


xcp = [];
if( iplot )
  subplot(1,2,2)
  hold on
end
for i_bdr = [5,6]
  ix = fea.grid.b(3,:) == i_bdr;
  ic = fea.grid.b(1,ix);
  ie = fea.grid.b(2,ix);
  je = mod(ie,size(fea.grid.c,1)) + 1;
  iv = fea.grid.c(sub2ind(size(fea.grid.c),ie,ic));
  jv = fea.grid.c(sub2ind(size(fea.grid.c),je,ic));
  p_bdr = fea.grid.p(:,[iv(1),jv]);
  l_choord = p_bdr(1,:)/cos(alfa*pi/180);

  col = 'b-';
  if( i_bdr==6 )
    col = 'r-';
  end
  cp  = evalexpr( 'cp', p_bdr, fea );
  xcp = [ xcp; [l_choord(:),cp(:)] ];
  if( iplot )
    plot( l_choord, cp, col, 'linewidth', 1.5 )
  end
end
if( iplot )
  axis([0 1 -5 1])
  set(gca,'Ydir','reverse')
  grid on
  ylabel('Pressure coefficient, c_p')
end

err = [];
if( strcmp(series,'0012') && alfa==6 )
  xycp_ref = l_xfoil_naca0012_alfa6();
  if( iplot )
    plot( xycp_ref(:,1), xycp_ref(:,3), 'k.' )
    legend( 'Lower surface', 'Upper surface', 'Xfoil' )
  end

  a     = polyarea( xcp(:,1), xcp(:,2) );
  a_ref = polyarea( xycp_ref(:,1), xycp_ref(:,3) );
  err = abs(a_ref-a)/a_ref;
elseif( iplot )
  legend( 'Lower surface', 'Upper surface' )
end

%------------------------------------------------------------------------------%
function [ xycp ] = l_xfoil_naca0012_alfa6()
% Pressure coeffcient cp computed with Xfoil for a NACA0012
% wing profile at 6 degrees angle of attack.

xycp = [ 1.00000  0.00126  0.42152;
         0.99372  0.00214  0.25098;
         0.98277  0.00365  0.18409;
         0.96999  0.00539  0.13218;
         0.95567  0.00731  0.09081;
         0.94026  0.00934  0.05619;
         0.92422  0.01142  0.02500;
         0.90784  0.01350 -0.00229;
         0.89128  0.01556 -0.02707;
         0.87462  0.01759 -0.05027;
         0.85793  0.01959 -0.07162;
         0.84120  0.02155 -0.09192;
         0.82446  0.02348 -0.11143;
         0.80771  0.02537 -0.13015;
         0.79095  0.02723 -0.14824;
         0.77419  0.02905 -0.16611;
         0.75742  0.03083 -0.18336;
         0.74066  0.03257 -0.20067;
         0.72389  0.03428 -0.21770;
         0.70712  0.03594 -0.23450;
         0.69035  0.03757 -0.25164;
         0.67359  0.03916 -0.26863;
         0.65683  0.04071 -0.28553;
         0.64008  0.04221 -0.30292;
         0.62333  0.04367 -0.32050;
         0.60659  0.04509 -0.33780;
         0.58987  0.04646 -0.35583;
         0.57315  0.04778 -0.37397;
         0.55645  0.04905 -0.39246;
         0.53976  0.05026 -0.41134;
         0.52309  0.05143 -0.43063;
         0.50644  0.05253 -0.45038;
         0.48980  0.05357 -0.47080;
         0.47319  0.05455 -0.49114;
         0.45660  0.05546 -0.51277;
         0.44004  0.05630 -0.53438;
         0.42351  0.05707 -0.55710;
         0.40701  0.05776 -0.58007;
         0.39055  0.05837 -0.60429;
         0.37412  0.05889 -0.62888;
         0.35774  0.05932 -0.65436;
         0.34140  0.05965 -0.68106;
         0.32511  0.05988 -0.70850;
         0.30888  0.06000 -0.73718;
         0.29272  0.06001 -0.76706;
         0.27662  0.05989 -0.79816;
         0.26060  0.05965 -0.83091;
         0.24467  0.05927 -0.86519;
         0.22884  0.05875 -0.90134;
         0.21313  0.05807 -0.93937;
         0.19755  0.05723 -0.98050;
         0.18214  0.05622 -1.02370;
         0.16693  0.05503 -1.07044;
         0.15196  0.05365 -1.12070;
         0.13730  0.05207 -1.17540;
         0.12303  0.05029 -1.23507;
         0.10929  0.04832 -1.30011;
         0.09622  0.04618 -1.37105;
         0.08399  0.04389 -1.44848;
         0.07278  0.04150 -1.53170;
         0.06271  0.03908 -1.62035;
         0.05380  0.03667 -1.71271;
         0.04603  0.03431 -1.80915;
         0.03929  0.03203 -1.90696;
         0.03346  0.02983 -2.00844;
         0.02840  0.02771 -2.10893;
         0.02399  0.02566 -2.21120;
         0.02013  0.02367 -2.31093;
         0.01674  0.02173 -2.40983;
         0.01376  0.01982 -2.50164;
         0.01114  0.01793 -2.58590;
         0.00883  0.01606 -2.65063;
         0.00682  0.01419 -2.69399;
         0.00508  0.01231 -2.69706;
         0.00360  0.01041 -2.64369;
         0.00238  0.00851 -2.52207;
         0.00141  0.00659 -2.30896;
         0.00070  0.00467 -2.00544;
         0.00025  0.00277 -1.62633;
         0.00003  0.00091 -1.20265;
         0.00003 -0.00091 -0.78171;
         0.00025 -0.00277 -0.36406;
         0.00070 -0.00467  0.02209;
         0.00141 -0.00659  0.35451;
         0.00238 -0.00851  0.61143;
         0.00360 -0.01041  0.79407;
         0.00508 -0.01231  0.91066;
         0.00682 -0.01419  0.97447;
         0.00883 -0.01606  0.99870;
         0.01114 -0.01793  0.99458;
         0.01376 -0.01982  0.97039;
         0.01674 -0.02173  0.93300;
         0.02013 -0.02367  0.88600;
         0.02399 -0.02566  0.83322;
         0.02840 -0.02771  0.77626;
         0.03346 -0.02983  0.71682;
         0.03929 -0.03203  0.65654;
         0.04603 -0.03431  0.59560;
         0.05380 -0.03667  0.53556;
         0.06271 -0.03908  0.47740;
         0.07278 -0.04150  0.42164;
         0.08399 -0.04389  0.37030;
         0.09622 -0.04618  0.32327;
         0.10929 -0.04832  0.28126;
         0.12303 -0.05029  0.24435;
         0.13730 -0.05207  0.21176;
         0.15196 -0.05365  0.18355;
         0.16693 -0.05503  0.15912;
         0.18214 -0.05622  0.13772;
         0.19755 -0.05723  0.11937;
         0.21313 -0.05807  0.10332;
         0.22884 -0.05875  0.08951;
         0.24467 -0.05927  0.07752;
         0.26060 -0.05965  0.06711;
         0.27662 -0.05989  0.05825;
         0.29272 -0.06001  0.05062;
         0.30888 -0.06000  0.04413;
         0.32511 -0.05988  0.03852;
         0.34140 -0.05965  0.03404;
         0.35774 -0.05932  0.03017;
         0.37412 -0.05889  0.02694;
         0.39055 -0.05837  0.02444;
         0.40701 -0.05776  0.02261;
         0.42351 -0.05707  0.02095;
         0.44004 -0.05630  0.02007;
         0.45660 -0.05546  0.01956;
         0.47319 -0.05455  0.01910;
         0.48980 -0.05357  0.01924;
         0.50644 -0.05253  0.01969;
         0.52309 -0.05143  0.02018;
         0.53976 -0.05026  0.02135;
         0.55645 -0.04905  0.02235;
         0.57315 -0.04778  0.02361;
         0.58987 -0.04646  0.02532;
         0.60659 -0.04509  0.02703;
         0.62333 -0.04367  0.02914;
         0.64008 -0.04221  0.03117;
         0.65683 -0.04071  0.03354;
         0.67359 -0.03916  0.03619;
         0.69035 -0.03757  0.03871;
         0.70712 -0.03594  0.04208;
         0.72389 -0.03428  0.04513;
         0.74066 -0.03257  0.04879;
         0.75742 -0.03083  0.05261;
         0.77419 -0.02905  0.05675;
         0.79095 -0.02723  0.06133;
         0.80771 -0.02537  0.06693;
         0.82446 -0.02348  0.07169;
         0.84120 -0.02155  0.07897;
         0.85793 -0.01959  0.08563;
         0.87462 -0.01759  0.09376;
         0.89128 -0.01556  0.10355;
         0.90784 -0.01350  0.11414;
         0.92422 -0.01142  0.12753;
         0.94026 -0.00934  0.14359;
         0.95567 -0.00731  0.16330;
         0.96999 -0.00539  0.18878;
         0.98277 -0.00365  0.22481;
         0.99372 -0.00214  0.27506;
         1.00000 -0.00126  0.42152];
