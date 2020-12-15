function [ c_paths ] = addpaths( i_add )
%ADDPATHS Add/remove app directories to/from working paths.
%
% [ C_PATHS ] = ADDPATHS( I_ADD ) Adds/removes app directories to/from
% working paths. The I_ADD (default true) indicates adding paths, and
% 0/false for removing paths. Optionally, returns a cell array with the
% added/removed paths in C_PATHS where the first column consists of
% all paths relative to the root app/mfile directory, and the second
% with corresponding absolute paths.

% Copyright 2013-2020 Precise Simulation, Ltd.

APP_NAME = 'featool';


fcn = @addpath;
if( nargin && ~i_add )
  fcn = @rmpath;
end

s_mfiledir = fileparts( which([mfilename('fullpath'),'.m']) );
s_paths = genpath( s_mfiledir );

fcn( s_paths );
s_app_path = l_get_app_path( APP_NAME );
fcn( s_app_path );

if( nargout )
  c_paths = regexp( strtrim(s_paths), ';+', 'split' );
  c_paths(cellfun(@isempty,c_paths)) = [];
  c_paths = [ regexprep(strrep(c_paths(:),s_mfiledir,''),['^',filesep()],''), c_paths(:) ];
  c_paths = [ c_paths; {'', s_app_path} ];
end


%------------------------------------------------------------------------------%
function [ s, s_home ] = l_get_app_path( app_name )
%L_GET_APP_PATH Duplicated function util/get_app_path

s_home = which( [app_name,'.p'] );
if( isempty(strtrim(s_home)) )
  s_home = which( [app_name,'.m'] );
end
s_home = fileparts(s_home);


s_userpath = strtrim( userpath() );
if( isempty(s_userpath) )
  userpath('reset');
  s_userpath = strtrim( userpath() );
end
while( length(s_userpath)>0 && ...
       ( s_userpath(end) == ';' || s_userpath(end) == ':' ) )
  s_userpath = strtrim(s_userpath(1:end-1));
end


if( isempty(s_userpath) )
  try
    add_folders = true;
    if( ispc() )
      s_userpath = winqueryreg( 'HKEY_CURRENT_USER',...
                                ['Software\Microsoft\Windows\CurrentVersion\' ...
                                   'Explorer\Shell Folders'], 'Personal' );
      s_userpath = strrep(s_userpath,[filesep,'Documents'],'');
    else
      s_userpath = char(java.lang.System.getProperty('user.home'));
    end
    if( isempty(strtrim(s_userpath)) )
      s_userpath = getenv('HOME');
    end
    if( isempty(strtrim(s_userpath)) )
      error( 'Could not find valid user path.' )
    end
  catch
    add_folders = false;
    s_userpath = s_home;
  end
  if( add_folders )
    s_userpath = fullfile(strtrim(s_userpath),'Documents','MATLAB');
  end
end

s = fullfile(strtrim(s_userpath),['.',lower(app_name)]);

if( exist(s)~=7 )
  st = rmkdir( s );
  if( st~=1 )
    warning( ['Could not find or create app home path in: ',char([10,10]),s] )
    s = '';
  end
end
