%EXPORT_VTP Export grid and data in VTP format.
%
%   [ S ] = IMPEXP_VTP( FILENAME, MODE, DATA, EXPR, SOLNUM )
%
%   Export of grid points and corresponding solution/expression
%   data in VTP ASCII format.
%
%   FILENAME specifies the output file name. MODE can only be a string
%   indicating export (default). DATA must be a valid grid or FEA data
%   struct. EXPR and SOLNUM are optional arguments for export, where
%   EXPR is a cell array with expressions to process. Expressions will
%   be evaluated in points corresponding to SCALARS/POINT_DATA format,
%   and dependent variables present in the FEA struct will be
%   prepended to EXPR. SOLNUM is optionally the solution number of
%   which to evaluate the expressions. FID_LOG is an optional log file
%   handle for message output (negative for gui output or empty for no
%   output).

% Copyright 2013-2023 Precise Simulation, Ltd.
