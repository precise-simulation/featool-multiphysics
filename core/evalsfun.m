%EVALSFUN Shape function evaluation driver subroutine.
%
%   [ VBASE, NLDOF, XLDOF, SFUN ] = EVALSFUN( SFUN, I_EVAL, N_SDIM, N_VERT, I_DOF, XI, AINVJAC, VBASE )
%   Evaluates local basis function I_DOF in point XI for shape function SFUN on all
%   cells corresponding to AINVJAC. Returns the evaluated values in VBASE and
%   the local number of degrees of freedom on vertices, edges, faces, and cell
%   interiors in NLDOF and XLDOF as well as the low level function name in SFUN.
%
%   See also EVALDVAR, EVALEXPR

% Copyright 2013-2022 Precise Simulation, Ltd.
