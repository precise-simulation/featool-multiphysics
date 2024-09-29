%TURBULENCE_INLETBCCALC Compute turbulence quantities for inlet boundaries.
%
%   [ KEO ] = TURBULENCE_INLETBCCALC( FEA, INTENSITY, LSCALE, I_MODE )
%   Computes turbulence quantities for inlet boundaries given a valid
%   FEA struct. INTENSITY specifies the turbulent intensity per inlet
%   boundary (default 0.1 = 10%), and LSCALE the fraction of the inlet
%   boundary length to compute the turbulence length scale (default
%   0.08). Returns a row array KEO with turbulent kinetic energy k,
%   and (specific) turbulent dissipation rate e/o for each inlet boundary.
%   I_MODE is the number of the physics mode to use.

% Copyright 2013-2024 Precise Simulation, Ltd.
