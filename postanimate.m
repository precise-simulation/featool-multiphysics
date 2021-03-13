function [ varargout ] = postanimate( varargin )
%POSTANIMATE Plot and animate solution.
%
%   [ M ] = POSTANIMATE( PROB, VARARGIN ) Postprocessing function to
%   plot and animate FEATool solutions. PROB is a valid finite element
%   problem struct. If an IMGNAME is given video or frames will be
%   saved in the corresponding IMGFMT (images will be appended with
%   the frame number before the image file extension). Optionally also
%   returns image frames in the struct M that can be played with the
%   movie command. Accepts the following property/value pairs as for
%   postplot and the following additional properties
%
%       Property    Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       start       integer/{1}               Solution to start animation
%       step        integer/{1}               Step size (positive/negative integer)
%       end         integer/{end}             Solution to end animation
%       pause       scalar/{0}                Pause between frame generation (seconds)
%       imgname     string/{}                 Image/video base filename
%       imgfmt      png,mp4,avi/{jpeg}        Output image/video format
%
% See also POSTPLOT, MOVIE, VIDEOWRITER, PRINT

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help postanimate, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'postanimate', varargin{:} );
if( ~nargout ), clear varargout; end
