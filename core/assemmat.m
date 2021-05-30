function [ varargout ] = assemmat( varargin )
%ASSEMMAT Assemble monolithic matrix.
%
%   [ A, T_A, T_SP ] = ASSEMMAT( PROB, S_A, I_CUB, N_CMAX, F_SPARSE, I_HRZ, SOLCOMP )
%   Called from assembleprob to assemble a monolithic matrix.
%
%   See also ASSEMBLEA

% Copyright 2013-2021 Precise Simulation, Ltd.

if( ~nargin && ~nargout ), help assemmat, return, end
varargout = cell( 1, nargout );
[varargout{:}] = featool( 'feval', 'assemmat', varargin{:} );
if( ~nargout ), clear varargout; end
