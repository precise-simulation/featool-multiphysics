%IMPORT_IMG Import and convert image to geometry.
%
%   [ GEOM ] = IMPORT_IMG( FILENAME, RES, SCALE, LEVEL, BORDER )
%   Import bitmap image with FILENAME image (JPG, PNG, BMP, TIFF),
%   extract contours, and convert to a geometry of closed polygons.
%
%   RES specifies the resolution in pixels of the polygon segments
%   (default 10 segments). SCALE is used for scaling from pixels to
%   geometry dimensions (default 1 = one geometry length unit per
%   pixel). LEVEL is used to prescribe an object detection threshold
%   using the image color/intensity level (0 < LEVEL < 1, default 0.5).
%   With the BORDER flag set an inverted border is added around the
%   image (unset by default).

% Copyright 2013-2022 Precise Simulation, Ltd.
