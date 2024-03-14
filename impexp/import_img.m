%IMPORT_IMG Import and convert image to geometry.
%
%   [ GEOM ] = IMPORT_IMG( FILENAME, TYPE, RES, SCALE, LEVEL, BORDER )
%   Import bitmap image with FILENAME image (JPG, PNG, BMP, TIFF),
%   extract contours, and convert to a geometry of closed polygons.
%
%   TYPE indicates reconstruction type (1 - polygon (default),
%   2 - spline curve, 3 - Bezier curve). RES specifies the resolution
%   in pixels of the curve segments (default 10 segments). SCALE is
%   used for scaling from pixels to geometry dimensions (default
%   1 - one geometry length unit per pixel). LEVEL is used to
%   prescribe an object detection threshold using the image
%   color/intensity level (0 < LEVEL < 1, default 0.5). With the
%   BORDER flag set an inverted border is added around the image
%   (unset by default).
%
%   Calling IMPORT_IMG without any input arguments opens the dialog
%   box with controls for tuning the image import.
%
%   Examples:
%
%      1) Open dialog box for image geometry import.
%
%      geom = import_img();
%
%      2) Convert image "filename" to geometry without dialog box.
%
%      geom = import_img( filename );

% Copyright 2013-2024 Precise Simulation, Ltd.
