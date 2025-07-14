%RUNTIMECLIENT MATLAB/GNU Octave TCP socket client.
%
%   RUNTIMECLIENT( VARARGIN ) Interpreter client to attach to a
%   socket. Reads, and evaluates received messages, and sends back
%   resulting output.
%
%   Interactive mode (default) allows use of the command line REPL
%   while the interpreter is pollilng for messages (every t_poll
%   seconds). In non-interactive mode reading (and writing) is blocked
%   until a STOP (127) terminated message has been received.
%
%   Accepts the following property/value pairs.
%
%       Property    Value        Description
%   -------------------------------------------------------------------
%       host        '127.0.0.1'  Socket server hostname or IPv4 address
%       port        3000         Socket server listening port
%       bufsize     65536        Socket buffer size
%       path        ''           Path to add to search paths
%       encoding    'utf-8'      Socket message string encoding
%       t_max       inf          Maximum time allowed for client
%       t_poll      0.05         Polling rate in interactive mode
%       t_wait      30.0         Maximum time allowed for startup
%       interactive true         Enable interactive mode and desktop
%       no_gui      false        Disable/hide desktop GUI
%       reconnect   3            Number of times to try to reconnect to client
%       quit        false        Quit interpreter after stopping/uncaught error
%       debug       0/{1,2}      Debugging message level (0 = off)

% Copyright 2013-2025 Precise Simulation, Ltd.











%------------------------------------------------------------------------------%
% Run client in RPC mode until a stop signal is received.











%------------------------------------------------------------------------------%
% Connect to server (blocking until t_wait).







%------------------------------------------------------------------------------%
% Find free local port.


%------------------------------------------------------------------------------%



% Receive request message.




% Evaluate request message.












% Send response message.


%------------------------------------------------------------------------------%





%------------------------------------------------------------------------------%



%------------------------------------------------------------------------------%


%------------------------------------------------------------------------------%



%------------------------------------------------------------------------------%


%------------------------------------------------------------------------------%


%------------------------------------------------------------------------------%


%------------------------------------------------------------------------------%









%------------------------------------------------------------------------------%




%------------------------------------------------------------------------------%


%------------------------------------------------------------------------------%








%------------------------------------------------------------------------------%




%------------------------------------------------------------------------------%
% Decode byte data and return message (version 1).

% Decode header.

% if( any(msg.mtype == [1,2]) )
%   msg.request = flags(4) == '1';
%   msg.response = flags(5) == '1';
% end

% Decode data.

%------------------------------------------------------------------------------%
% Encode message and return byte data.


% Encode data.

% Encode header.
% flags = '00100000';
% if( request ), flags(4) = '1'; end
% if( response ), flags(5) = '1'; end


%------------------------------------------------------------------------------%









%------------------------------------------------------------------------------%











%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%




%------------------------------------------------------------------------------%





%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%




