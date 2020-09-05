function [L, U, levs] = iluk(A, lfil, Lvl)
%ILUK Incomplete LU-factorization with fill level k
%
%       [L, U] = iluk(A, lfil) computes an Incomplete LU factorization of the
%   matrix A with level of fill lfil. The matrix A is expected to be
%   sparse, real, and square. The parameter lfil can be any nonnegative
%   integer. L and U are the lower and upper triangular matrices,
%   respectively, of the incomplete factorization.
%
%   The algorithm proceeds in two phases. In the first phase (symbolic)
%   locations of permitted fill-in entries are determined based on the
%   level of fill parameter lfil and are stored in a data structure. In the
%   second phase (numeric) the numeric values of permitted fill-in entries
%   are computed and the L and U factors are constructed.
%
%   [L, U, levs] = iluk(A, lfil) additionally returns sparse matrix levs,
%   which is the TRANSPOSED level of fill matrix. That is, levs(i,j) is the
%   level of fill for the (j,i) element of A. Note that the levels of fill
%   in levs are shifted by 1 so that any entries with level of fill 1
%   correspond to nonzero entries in A, i.e., level 0 entries.
%
%   If the iluk method is being repeatedly called such as in a loop and if
%   the sparsity pattern of A is not changing between calls then the
%   following can be done to save computations. In the first iteration call
%   [L, U, levs] = iluk(A, lfil) and keep the levs matrix. In subsequent
%   iterations call [L, U] = iluk(A, [], levs), which avoids doing a
%   symbolic factorization by using levs to identify permitted fill-in
%   entries. Such usage is typical, for example, in time-dependent PDE
%   problems.
%
%   Before using iluk you must select a C++ compiler in MATLAB using
%   mex -setup. Once the compiler has been selected you need to issue
%   the following commands to compile the necessary files:
%
%   mex -largeArrayDims iluk_Cmex.cpp iluk_helper.cpp
%   mex -largeArrayDims iluk_Cmex_nosymb.cpp iluk_helper.cpp

% Copyright (c) 2014 Killian Miller. All rights reserved.

if(nargin < 2 || nargin > 3)
  error('Wrong number of inputs');
end

if(isempty(A) || ~isreal(A) || ~issparse(A) || size(A,1) ~= size(A,2))
  error('Matrix A must be nonempty, square, sparse, and real');
end

if(nargin == 2 && (lfil < 0 || lfil ~= floor(lfil)))
  error('Level of fill must be a nonnegative integer');
end

if(nargin == 3 && nargout == 3)
  error('Only two outputs permitted with three inputs');
end

if(nargin == 3)
  if(isempty(Lvl) || ~isreal(Lvl) || ~issparse(Lvl) || size(Lvl,1) ~= size(Lvl,2))
    error('Matrix Lvl must be nonempty, square, sparse, and real');
  end
  if(size(Lvl,1) ~= size(A,1))
    error('Matrix Lvl must be the same size as A');
  end
end

if(nargout < 2)
  error('At least two outputs required');
elseif(nargout == 2)
  A = A';
  if(nargin == 2)
    [L, U] = iluk_Cmex(A, lfil);
  elseif(nargin == 3)
    [L, U] = iluk_Cmex_nosymb(A, Lvl);
  end
elseif(nargout == 3)
  A = A';
  [L, U, levs] = iluk_Cmex(A, lfil);
else
  error('Too many outputs');
end
