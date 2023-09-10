%ASSEMMAT Assemble monolithic matrix.
%
%   [ A, T_A, T_SP ] = ASSEMMAT( PROB, S_A, I_CUB, N_CMAX, F_SPARSE, I_HRZ, SOLCOMP )
%
%   Assemble monolithic (coupled) matrix (for all dependent variables)
%   for the field S_A in the finite element data struct PROB. S_A can
%   be either "m" or "a", designating assembly of the mass matrix
%   defined in PROB.EQN.M, or system matrix defined in PROB.EQN.A.
%
%       Input       Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       prob        struct                 FEA problem data struct
%       s_a         char    {m/a}          Matrix field to assemble (mass/system)
%       i_cub       scalar  {2}            Numerical integration rule
%       n_cmax      scalar  {50000}        Max number of cells to assemble
%                                          for at once (to limit memory consumption)
%       f_sparse    logical {true}         Return sparse/struct matrix format
%       i_hrz       logical {false}        Apply HRZ diagonal (mass) lumping
%       solcomp     {all dvars/subd}       Dependent variables/subdomains to assemble for
%                                                                                         .
%       Output      Value/[Size]           Description
%       -----------------------------------------------------------------------------------
%       A           sparse/struct [n_A]    Assembled system matrix
%       t_a         scalar                 Time spent assembling matrix
%       t_sp        scalar                 Time for sparse matrix conversion
%
%   See also ASSEMBLEA, ASSEMBLEPROB

% Copyright 2013-2023 Precise Simulation, Ltd.
