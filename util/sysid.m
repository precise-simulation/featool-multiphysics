function [ id ] = sysid( type )
%SYSID Return system id (12 [0-9A-F] chars).

% Copyright 2013-2021 Precise Simulation, Ltd.
if( ~nargin ), type = 2; end

id = '';

% Sysid type 1.
if( type < 2 )
  try
    [~,key_ver] = get_license_key();
    if( key_ver == 1 )
      id = l_getmac1();
      return;
    end
  catch,end
end

% Sysid type 2.
if( ispc() )

  try
    [st,id] = builtin( 'system', 'wmic path win32_computersystemproduct get uuid' );
    if( ~isequal(st,0) ), error( 'Can not get Windows system id.' ), end

    id(id == ':' | id == '-' | isspace(id)) = ''; id = upper(id(end-11:end));
  catch, id = ''; end
  if( all(id(1)==id) ), id = ''; end

elseif( ismac() )

  try
    [st,id] = builtin( 'system', 'system_profiler SPHardwareDataType' );
    if( ~isequal(st,0) ), error( 'Can not get Mac system id.' ), end

    id = regexp(upper(id), 'HARDWARE UUID:\s*([0-9A-F|\-].*)', 'tokens', 'once');
    if( iscell(id) ), id = id{1}; end
    id(id == ':' | id == '-' | isspace(id)) = ''; id = id(end-11:end);
  catch, id = ''; end

end

if( ~isequal(1:12, regexp(id,'[0-9A-F]')) )
  try
    id = l_get_mac();
  catch, id = ''; end
end

if( ~(isequal(1:12, regexp(id,'[0-9A-F]')) || isequal(1:17, regexp(id,'[0-9A-F]'))) || ...
    all(id == '0') || all(id == 'F') )
  % Assign random id with end/error char set to Z.
  id = dec2hex(randi(15,1,12))';
  id(end) = 'Z';
end
id = id(1:12);


%------------------------------------------------------------------------------%
function [ mac ] = l_get_mac()

try

  mac = {};
  ni = java.net.NetworkInterface.getNetworkInterfaces;
  while( ni.hasMoreElements )
    addr = ni.nextElement.getHardwareAddress;
    if( ~isempty(addr) )
      mac = [ mac, {sprintf('%.2X', typecast(addr,'uint8'))} ];
    end
  end
  for i=1:length(mac)
    if( ~isequal(1:12, regexp(upper(mac{i}),'[0-9A-F]')) )
      mac{i} = '';
    end
  end

catch

  if( ispc() )
    macfunc = 'getmac';
    % elseif( isunix() )
    %   if( system('grep -q Microsoft /proc/version')==0 )   % For WSL bash.
    %     macfunc = '/mnt/c/Windows/System32/getmac.exe';
    %   end
  else
    macfunc = 'ifconfig';
    [st,res] = system(['which ',macfunc]);
    if( st~=0 || isempty(strfind(res,macfunc)) )
      macfunc = 'ip a';
    end
  end

  [st,mac] = system( macfunc );
  if( isequal(st,0) )
    mac = l_sysid_str2mac( mac );
  end

end

if( ~iscell(mac) ), mac = {mac}; end
mac = unique(mac);
for i=1:length(mac)
  mac_i = mac{i}; if( iscell(mac_i) ), mac_i = mac_i{1}; end
  if( ~(all(mac_i == '0') || all(mac_i == 'F')) )
    mac = mac_i;
    break;
  elseif( i == length(mac) )
    mac = '';   % No mac found.
  end
end


%------------------------------------------------------------------------------%
function [ mac ] = l_getmac1()

try
  if( ispc() )
    macfunc = 'getmac';
    % elseif( isunix() )
    %   if( system('grep -q Microsoft /proc/version')==0 )   % For WSL bash.
    %     macfunc = '/mnt/c/Windows/System32/getmac.exe';
    %   end
  else
    macfunc = 'ifconfig';
    [st,res] = system(['which ',macfunc]);
    if( st~=0 || isempty(strfind(res,macfunc)) )
      macfunc = 'ip a';
    end
  end

  try
    [tmp,mac] = builtin( 'system', macfunc );
  catch
    [tmp,mac] = system( macfunc );
  end

  mac = l_sysid_str2mac( mac );

  mac = unique(mac);
  for i=1:length(mac)
    mac_i = mac{i}; if( iscell(mac_i) ), mac_i = mac_i{1}; end
    if( i==length(mac) || ~(all(mac_i == '0') || all(mac_i == 'F')) )
      mac = mac_i;
      break;
    end
  end

  if( ~ischar(mac) || isempty(mac) || ~(length(mac)==12 || length(mac)==17) )
    % fprintf( fid, '%s\n', [datestr(now()),': ERROR : sysinfo-net1 : ',mac,' failed validation'] );
    % fprintf( fid, '%s\n', [datestr(now()),': ',lasterr()] );
    mac = dec2hex(randi(15,1,12))';
    mac(end) = 'Z';
    % ier = 1;   % Mac verification failure.
  end

catch
  % fprintf( fid, '%s\n', [datestr(now()),': ERROR : sysinfo-net2'] );
  % fprintf( fid, '%s\n', [datestr(now()),': ',lasterr()] );
  mac = dec2hex(randi(15,1,12))';
  mac(end) = 'Z';
  % ier = 2;   % Mac syscall failure.
end

%------------------------------------------------------------------------------%
function [ cmac ] = l_sysid_str2mac( mac )

mac = regexp( mac, '\n|\r', 'split' );
ind_docker = find(~cellfun(@isempty, regexp(lower(mac), 'docker\d:')));

mac = regexp( upper(mac), '([0-9A-F]{2}[-|:]){5}[0-9A-F]{2}', 'match' );
ind_mac = find(~cellfun(@isempty, mac));

for i=1:length(ind_docker)
  i_del = ind_mac(find(ind_mac > ind_docker(i), 1));
  mac{i_del} = [];
end
mac = mac(~cellfun(@isempty,mac));

cmac = {};
for i=1:length(mac)
  mac_i = mac{i};
  if( iscell(mac_i) )
    cmac = [ cmac, mac_i(:)' ];
  else
    cmac = [ cmac, {mac_i} ];
  end
end
mac = cmac;

for i=1:length(cmac)
  mac_i = cmac{i}; if( iscell(mac_i) ), mac_i = mac_i{1}; end
  mac_i(mac_i == ':' | mac_i == '-' | isspace(mac_i)) = '';
  cmac{i} = mac_i;
end
