function [ fea, out ] = ex_classicpde2( varargin )
%EX_CLASSICPDE2 Eigenmodes for a L-shaped membrane.
%
%   [ FEA, OUT ] = EX_CLASSICPDE2( VARARGIN ) Eigenmodes for a
%   L-shaped membrane. Accepts the following property/value pairs.
%
%       Input       Value/{Default}        Description
%       -----------------------------------------------------------------------------------
%       hmax        scalar {0.05}          Grid cell size
%       sfun        string {sflag1}        Shape function
%       iplot       scalar 0/{1}           Plot solution (=1)
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.


cOptDef = { ...
  'hmax',     0.05; ...
  'sfun',     'sflag1'; ...
  'iplot',    1; ...
  'tol',      0.05; ...
  'fid',      1 };
[got,opt] = parseopt(cOptDef,varargin{:});
fid       = opt.fid;


% Geometry definition.
gobj = gobj_polygon( [-1  1 1 0 0 -1;
                      -1 -1 1 1 0  0].', 'P1' );
fea.geom.objects = { gobj };


% Grid generation.
fea.grid = gridgen(fea,'hmax',opt.hmax,'fid',fid);
n_bdr = max(fea.grid.b(3,:));   % Number of boundaries.


% Problem definition.
fea.sdim  = { 'x' 'y' };   % Coordinate names.

fea = addphys( fea, @poisson, {'u'} );
fea.phys.poi.sfun  = { opt.sfun };

fea.phys.poi.bdr.coef{1,end}  = repmat({0},1,n_bdr);

fea = parsephys(fea);


% Parse and solve problem.
fea = parseprob(fea);
[fea.sol.u,fea.sol.l] = solveeig( fea, 'fid', fid, 'neigs', 20 );


% Postprocessing.
if( opt.iplot>0 )
  postplot( fea, 'surfexpr', 'u', 'surfhexpr', 'u', 'solnum', 12 )
  title(['Solution lambda 12 = ',num2str(fea.sol.l(12))])
end


l_ref = [9.640368;15.197253;19.739209;29.521482;31.914209;41.475693;44.948492;49.348023;49.348023;56.710931;65.376542;71.059312;71.572682;78.956838;89.306363;92.306914;97.380734;98.69605;101.607501;112.369198];
out.err  = norm( l_ref - fea.sol.l )/norm(l_ref);
out.pass = out.err < opt.tol;

if( ~nargout )
  clear fea out
end
