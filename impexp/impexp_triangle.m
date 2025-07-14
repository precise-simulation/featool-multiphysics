%IMPEXP_TRIANGLE Import/export grid data in Triangle format.
%
%   [ DATA ] = IMPEXP_TRIANGLE( FILE_NAME, MODE, DATA, HMAX, HMAX_E, FID_LOG ) Import
%   or export of Triangle grid data in .node/.ele/.edge/.poly format. FILE_NAME is a string
%   specifying the (root) file name to process. MODE can either be a string indicating
%   import or export. For export, DATA can be either a whole fea struct or just the geom
%   struct which is used to compute the poly information used by Triangle. HMAX and HMAX_E
%   specifies grid sizes for subdomains and edges when exporting. A DATA grid struct is
%   output when importng. FID_LOG is an optional log file for message output.

% Copyright 2013-2025 Precise Simulation, Ltd.
