%OPENFOAM_WRITE Write OpenFOAM dictionary.
%
%   OPENFOAM_WRITE( FILE_PATH, DICT, FID_LOG ) Write OpenFOAM dictionary.
%
%   DICT is a OpenFOAM dictionary struct with (nested) fields:
%
%       FoamFile (mandatory dictionary field)
%           version 2.0
%           format  ascii
%           class   dictionary
%           object  name
%       keyword1 argument1
%              ...
%       keywordN argumentN
%
%   which is written to the (file) FILE_PATH in ASCII format.
%
%   FID_LOG is a optional file identifier for terminal and message
%   output 1>=ouput to file (default 1/stdout), -1=output to gui
%   terminal, []=no output.
%
%   Examples:
%
%      1) Example of writing a transportModel dictionary.
%
%      dict.transportModel.FoamFile = struct('version', '2.0', 'format', 'ascii', 'class', 'dictionary', 'object', 'transportProperties');
%      dict.transportModel.transportModel = 'Newtonian';
%      dict.transportModel.nu = '[0 2 -1 0 0 0 0] 0.001';
%      casedir = 'C:\temp\oftest';
%      openfoam_write( casedir, dict )
%
%   See also OPENFOAM_READ

% Copyright 2013-2026 Precise Simulation, Ltd.
