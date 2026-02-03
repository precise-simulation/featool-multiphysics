%EXPORT_CSV Export grid and data in CSV format.
%
%   [ DATA ] = EXPORT_CSV( FILENAME, EXPR, XP, PROB, SOLNUM )
%
%   Evaluate expressions (EXPR) and export data in CSV format. The CSV
%   file will consist of a header line with the space coordinate names
%   (P), and given expressions (as columns). This is followed by the
%   actual data rows where one line corresponds to the data for the
%   specific grid point.
%
%   FILENAME specifies the output file name. EXPR is a string or cell
%   array with expressions to evaluate in coordinates XP (grid
%   coordinates if empty). PROB must be a valid FEA problem data
%   struct, SOLNUM is optionally the solution number of which to
%   evaluate the expressions.
%
%   Returns the evaulated DATA (with XP appended first).
%
%   See also EVALEXPR, EVALEXPRP

% Copyright 2013-2026 Precise Simulation, Ltd.
