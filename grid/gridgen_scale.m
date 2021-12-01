function [ varargout ] = gridgen_scale( varargin )
%GRIDGEN_SCALE Grid generation with scaling.
%
%   [ GRID, STATS ] = GRIDGEN_SCALE( SIN, VARARGIN ) Generates a
%   scaled grid/mesh for the geometry defined by the objects in SIN by
%   calling an external grid generation algorithm. SIN can be either a
%   valid fea problem struct with geometry defined in
%   SIN.geom.objects, a cell array of multiple geometry objects, or a
%   single geometry object. The geometry is first scaled according to
%   the input arguments, the scaled geometry is then meshed, and
%   finally the mesh the re-scaled to fit the original
%   geometry. Accepts the following property/value pairs are available
%   with the default grid generation algorithm.
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       axis        int {1}                   Scaling axis (1=x, 2=y, 3=z)
%       p_scal      scalar/array {midpoint}   Scaling points (along axis)
%       s_scal      scalar/array {[]}         Scaling factors (per scaling point)
%       varargin    pv pairs {defaults}       Property/value pairs to pass to gridgen
%
%   P_SCAL are coordinates along the AXIS to identify geometry objects
%   to scale. Geometry objects for which P_SCAL falls within the
%   bounding boxes will be scaled with the corresponding S_SCAL entry.
%
%   Examples:
%
%      1) 3x1 rectangle scaled by 1/3 in the x-direction.
%
%      grid = gridgen_scale( gobj_rectangle(0,3), 's_scal', 1/3 );
%      plotgrid( grid )
%
%      2) Polygon scaled by 1/5 in the y-direction.
%
%      grid = gridgen_scale( gobj_polygon([0 0;1 0;0.8 2;1.1 5; 0 5]), 'axis', 2, 's_scal', 1/5 );
%      plotgrid( grid )
%
%      3) Cylinder scaled by 1/6 in the z-direction.
%
%      grid = gridgen_scale( gobj_cylinder([0,0,0],0.5,5,3), 'axis', 3, 's_scal', 1/6 );
%      plotgrid( grid )
%
%      4) Scaling three parallel cylinders in a block (only scaling section with cylinders).
%
%      b1 = gobj_block(-0.75,3.25,-0.75,0.75,-0.5,5.5);
%      c1 = gobj_cylinder([0,0,0],0.5,5,3);
%      c2 = gobj_cylinder([1.25,0,0],0.5,5,3);
%      c3 = gobj_cylinder([2.5,0,0],0.5,5,3);
%      geom.objects = {b1, c1, c2, c3};
%      grid = gridgen_scale( geom, 'axis', 3, 's_scal', 1/5 );
%      plotgrid( grid, 'selcells', 'y>0' )
%      plotsubd( grid, 'linewidth', 1, 'labels', 'off' )
%      plotedg( grid, 'colors', {'k'}, 'labels', 'off', 'linewidth', 1 );
%      view(3)
%
%      4) Scale two blocks in a cylinder (only scaling sections with blocks).
%
%      s1 = gobj_sphere();
%      b1 = gobj_block(-1,1,-1,1,1.5,4.5);
%      b2 = gobj_block(-1,1,-1,1,-3.5,-1.5,'B2');
%      c1 = gobj_cylinder([0,0,-4],1.75,9,3);
%      geom.objects = {s1, b1, b2, c1};
%      grid = gridgen_scale( geom, 'axis', 3, 'p_scal', [-2.5, 3], 's_scal', [1/2, 1/4] );
%      plotgrid( grid, 'selcells', 'y>0' )
%      plotsubd( grid, 'linewidth', 1, 'labels', 'off' )
%      plotedg( grid, 'colors', {'k'}, 'labels', 'off', 'linewidth', 1 );
%      view(3)
%
%   See also GRIDGEN

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help gridgen_scale, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'gridgen_scale', varargin{:} );
if( ~nargout ), clear varargout; end
