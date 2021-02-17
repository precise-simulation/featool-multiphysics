function [ fea, out ] = ex_classicpde1( varargin )
%EX_CLASSICPDE1 Eigenmodes for a circular drum.
%
%   [ FEA, OUT ] = EX_CLASSICPDE1( VARARGIN ) Eigenmodes for a
%   circular drum. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       igrid       scalar 0/{1}           Cell type (0=quadrilaterals, 1=triangles)
%       hmax        scalar {0.1}           Grid cell size
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'igrid',    1; ...
  'hmax',     0.1; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      0.02; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_circle();
fea.geom.objects = { gobj };


% Grid generation.
if( opt.igrid==1 )
  fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
else
  fea.grid = circgrid( 16, 12, 1 );
  if( opt.igrid<0 )
    fea.grid = quad2tri( fea.grid );
  end
end
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };   % Coordinate names.

fea = addphys( fea, @poisson, {'u'} );
fea.phys.poi.sfun  = { opt.sfun };

fea.phys.poi.bdr.coef{1,end}  = repmat({0},1,n_bdr);

fea = parsephys(fea);


% Parse and solve problem.
fea = parseprob(fea);
[fea.sol.u,fea.sol.l] = solveeig( fea, 'fid', fid );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'u', 'surfhexpr', 'u' )
  title(['Solution lambda 6 = ',num2str(fea.sol.l(end))])
end


l_ref = [5.783186;14.681971;14.681971;26.374616;26.374617;30.471262];
f = sqrt(fea.sol.l)/(2*pi);
out.err  = [ norm( l_ref - fea.sol.l )/norm(l_ref) ;
             norm(abs(f([2,4,6])/f(1)-[1.59;2.14;2.30])./[1.59;2.14;2.30]) ];
out.pass = all(out.err < opt.tol);

if( ~nargout )
  clear fea out
end
