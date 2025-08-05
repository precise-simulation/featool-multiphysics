%GRIDFACE Generate grid cell face numbering.
%
%   [ F, FV, IEXT ] = GRIDFACE( C, B ) Generates an array F with global face
%   numbers and FV with grid point indices per face. Input data is an
%   array C with cell connectivities, and optionally the boundary
%   struct B if upper triangular (OpenFOAM) ordering should be used
%   (default normal ordering).
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       c           (n_v,n_c)              Cell connectivity pointer array where the
%                                          entries in each column correspond to the
%                                          cell vertex numbers specifying the cells.
%       b           struct                 Boundary struct to use upper triangular ordering
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (n_lf,n_c)             Array with global face numbers for all cells.
%       fv          (n_f,3/4)              Array with indices to vertices for each face.
%       iext        (n_f,1)                Logical array designating external faces .

% Copyright 2013-2025 Precise Simulation, Ltd.
