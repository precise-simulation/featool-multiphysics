%DECODE_JSON Decode JSON string to data struct.
%
%   [ DATA, IER ] = DECODE_JSON( JSON, USE_BUILTIN ) Parses a JSON
%   string and returns a struct with the parsed data. JSON objects are
%   converted to structures (USE_BUILTIN default true uses jsondecode
%   if available).
%
%   Example:
%
%     google_search  = 'http://ajax.googleapis.com/ajax/services/search/web?v=1.0&q=matlab';
%     matlab_results = decode_json( urlread(google_search) );
%     disp( matlab_results{1}.responseData.results{1}.titleNoFormatting )
%     disp( matlab_results{1}.responseData.results{1}.visibleUrl )

% Based on https://www.mathworks.com/matlabcentral/fileexchange/20565-json-parser?focused=98def053-f9e3-da18-f7f1-7a862d818fbc with BSD Type License.

% Copyright 2013-2019 Precise Simulation Ltd.









%------------------------------------------------------------------------------%










%------------------------------------------------------------------------------%












%------------------------------------------------------------------------------%







%------------------------------------------------------------------------------%




%------------------------------------------------------------------------------%






%------------------------------------------------------------------------------%




%% UNIT TESTS

%!assert( isequal( decode_json('[]'), [] ) )
%!assert( isequal( decode_json('0.1'), 0.1 ) )
%!assert( isequal( decode_json('[1,2]'), [1;2] ) )
%!assert( isequal( decode_json('[[2,1]]'), [2,1] ) )
%!assert( isequal( decode_json('[[1,2],[2,3],[3,4]]'), [1,2;2,3;3,4] ) )
%!assert( isequal( decode_json('[[1,2,3],[4,5,6]]'), [1,2,3;4,5,6] ) )
%!assert( isequal( decode_json('["a","b"]'), {'a';'b'} ) )
%!assert( isequal( decode_json('[[[0,0],[1,0],[0,2]],[4,5]]'), {[0,0;1,0;0,2];[4;5]} ))
%!assert( isequal( decode_json('[[2,1],"a"]'), {[2;1];'a'} ) )
%!assert( isequal( decode_json('[[1,2],["c","d"],[3,"f"]]'), {[1;2];{'c';'d'};{3;'f'}} ) )
%!assert( isequal( decode_json('{"a":"b"}'), struct('a','b') ) )
%!assert( isequal( decode_json('["a",2,[]]'), {'a';2;[]} ) )
%!assert( isequal( decode_json('[["a",2,[]],[]]'), {{'a';2;[]};[]} ) )
