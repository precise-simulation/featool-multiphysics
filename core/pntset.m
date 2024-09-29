%PNTSET Set Dirichlet point constraints.
%
%   [ F, INDROW, AMAT, T ] = PNTSET( PROB, F, AMAT, ISYMM, SET_NULL )
%   Sets Dirichlet point constraints conditions in the right hand side
%   load vector F and global matrix AMAT with the information in the finite
%   element problem struct PROB.
%
%       Input       Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 Finite element problem struct
%       f           (neq,1)                Right hand side/load vector
%       amat        (n_a,n_a)              System matrix (sparse or triplet format)
%       isymm       scalar/{0}             Symmetrize BCs if applicable.
%       set_null    scalar/{0}             Set zeros in f vector.
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       f           (neq,1)                Modified right hand side/load vector
%       indrow      (neq,1)                Index to rows (dofs) in rhs which were set
%       amat        (n_a,n_a)              Modified system matrix
%       t           scalar                 Time spent in function

% Copyright 2013-2024 Precise Simulation, Ltd.
