function [ ux_x ] = ex_poisson3_derivative( x );
%EX_POISSON3_DERIVATIVE Derivative interpolation function
%
%   [ UX_X ] = EX_POISSON3_DERIVATIVE( X ) Derivative interpolation function

% Copyright 2013-2025 Precise Simulation, Ltd.

ux   = [0 0.107535709861634 0.170943549921419 0.217687143600172 0.253587789782147 0.281383488256056 0.302649539007050 0.318378206032763 0.329213990472178 0.335564675720330 0.337657240836525 0.335564675720314 0.329213990472148 0.318378206032723 0.302649539007007 0.281383488256008 0.253587789782100 0.217687143600132 0.170943549921393 0.107535709861621 0];
ux_x = interp1( linspace(0,1,numel(ux)), ux, x, 'pchip' );
ux_x( x<=0 | x>=1 ) = 0;
