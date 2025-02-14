%ADDPHYS Parse and add physics mode to problem struct.
%
%   [ PROB ] = ADDPHYS( PROB, F_PHYS, C_DVAR, F_OUT ) Parses and adds
%   the physics mode with function name or handle F_PHYS to the
%   problem struct PROB. C_DVAR prescribes the given instead of
%   default dependent variables names, and the F_OUT flag adds the
%   mode to the PROB struct (default true) or optionally simply
%   returns the computed physics mode struct.
%
%   Example:
%
%      1) Add default Poisson physics mode.
%
%      fea.sdim = {'x', 'y'};
%      fea.geom.objects = {gobj_circle};
%      fea = addphys( fea, @poisson );
%
%      2) Add convection and diffusion physics mode with 3 named dependent variables.
%
%      fea.sdim = {'x', 'y', 'z'};
%      fea.geom.objects = {gobj_block};
%      fea = addphys( fea, @convectiondiffusion, {'c1', 'c2', 'c3'} );
%
%   See also PARSEPHYS

% Copyright 2013-2025 Precise Simulation, Ltd.
