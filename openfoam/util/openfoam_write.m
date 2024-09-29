%OPENFOAM_WRITE Write OpenFOAM dictionary.
%
%   OPENFOAM_WRITE( PATH_, DICT, FID_LOG ) Write OpenFOAM dictionary.
%
%   DICT is a OpenFOAM dictionary struct with (nested) fields:
%
%       class (dictionary or list)
%       keyword (dictionary keyword or list name/number of items)
%       data (entry data)
%
%   which is written to the (file) PATH.
%
%   FID_LOG is a optional file identifier for terminal and message
%   output 1>=ouput to file (default 1/stdout), -1=output to gui
%   terminal, []=no output.
%
%   See also OPENFOAM_READ

% Copyright 2013-2024 Precise Simulation, Ltd.
