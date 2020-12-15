function [id]=dmumps(id,mat)
%
% [id]=dmumps(id,mat)
% id is a structure (see details in initmumps.m and MUMPS documentation)
% mat is optional if the job is -1 or -2
% mat is a square sparse matrice
% information are return in id fields
%
%****************************************
%This help menu gives details about the use of dmumps and initmumps
%****************************************
%
%--------------- Input Parameters ---------------
%
% - mat: sparse matrix which has to be provided as the second argument of dmumps if id.JOB is strictly larger than 0.
%
% - id.SYM: controls the matrix type (symmetric positive definite, symmetric indefinite or unsymmetric) and it has do be initialized by the user before the initialization phase of MUMPS (see id.JOB). Its value is set to 0 after the call of initmumps.
%
% - id.JOB: defines the action that will be realized by MUMPS: initialize, analyze and/or factorize and/or solve and release MUMPS internal C/Fortran data. It has to be set by the user before any call to MUMPS (except after a call to initmumps, which sets its value to -1).
%
% - id.ICNTL and id.CNTL: define control parameters that can be set after the initialization call (id.JOB = -1). See Section ``Control parameters'' of the MUMPS user's guide for more details. If the user does not modify an entry in id.ICNTL then MUMPS uses the default parameter. For example, if the user wants to use the AMD ordering, he/she should set id.ICNTL(7) = 0. Note that the following parameters are inhibited because they are automatically set within the interface: id.ICNTL(19) which controls the Schur complement option and id.ICNTL(20) which controls the format of the right-hand side. Note that parameters id.ICNTL(1:4) may not work properly depending on your compiler and your environment. In case of problem, we recommand to swith printing off by setting id.ICNL(1:4)=-1.
%
% - id.PERM\_IN: corresponds to the given ordering option (see Section ``Input and output parameters'' of the MUMPS user's guide for more details). Note that this permutation is only accessed if the parameter id.ICNTL(7) is set to 1.
%
% - id.COLSCA and id.ROWSCA: are optional scaling arrays (see Section ``Input and output parameters'' of the MUMPS user's guide for more details)
%
% - id.RHS: defines the right-hand side. The parameter id.ICNTL(20) related to its format (sparse or dense) is automatically set within the interface. Note that id.RHS is not modified (as in MUMPS), the solution is returned in id.SOL.
%
% - id.VAR\_SCHUR: corresponds to the list of variables that appear in the Schur complement matrix (see Section ``Input and output parameters'' of the MUMPS user's guide for more details).
%
% - id.REDRHS(input parameter only if id.VAR\_SCHUR was provided during the factorization and if ICNTL(26)=2 on entry to the solve phase): partial solution on the variables corresponding to the Schur complement. It is provided by the user and normally results from both the Schur complement and the reduced right-hand side that were returned by MUMPS in a previous call. When ICNTL(26)=2, MUMPS uses this information to build the solution id.SOL on the complete problem. See Section ``Schur complement'' of the MUMPS user's guide for more details.
%
%--------------- Output Parameters ---------------
%
% - id.SCHUR: if id.VAR\_SCHUR is provided of size SIZE\_SCHUR, then id.SCHUR corresponds to a dense array of size (SIZE\_SCHUR,SIZE\_SCHUR) that holds the Schur complement matrix (see Section ``Input and output parameters'' of the MUMPS user's guide for more details). The user does not have to initialize it.
%
% - id.REDRHS(output parameter only if ICNTL(26)=1 and id.VAR\_SCHUR was defined): Reduced right-hand side (or condensed right-hand side on the variables associated to the Schur complement). It is computed by MUMPS during the solve stage if ICNTL(26)=1. It can then be used outside MUMPS, together with the Schur complement, to build a solution on the interface. See Section ``Schur complement'' of the MUMPS user's guide for more details.
%
% - id.INFOG and id.RINFOG: information parameters (see Section ``Information parameters'' of the MUMPS user's guide ).
%
% - id.SYM\_PERM: corresponds to a symmetric permutation of the variables (see discussion regarding ICNTL(7) in Section ``Control parameters'' of the MUMPS user's guide ). This permutation is computed during the analysis and is followed by the numerical factorization except when numerical pivoting occurs.
%
% - id.UNS\_PERM: column permutation (if any) on exit from the analysis phase of MUMPS (see discussion regarding ICNTL(6) in Section ``Control parameters'' of the MUMPS user's guide ).
%
% - id.SOL: dense vector or matrix containing the solution after MUMPS solution phase. Also contains the nullspace in case of null space computation, or entries of the inverse, in case of computation of inverse entries.
%
%--------------- Internal Parameters ---------------
%
% - id.INST: (MUMPS reserved component) MUMPS internal parameter.
%
% - id.TYPE: (MUMPS reserved component) defines the arithmetic (complex or double precision).
%

errmsg = nargoutchk(1,1,nargout);
if(~isempty(errmsg))
     disp(errmsg);
     return;
end

arithtype = 1;

if(id.JOB == -2)
     if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
     end
     if(id.TYPE ~= arithtype)
       disp('You are trying to call z/d version on a d/z instance');
       return;
     end
     dmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS,id.KEEP,id.DKEEP);
     id = [];
     return;
end


if(id.JOB == -1)
     if(id.INST~=-9999)
         disp('Already initialized instance');
         return;
     end
     [inform,rinform,sol,inst,schur,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl,colsca_out,rowsca_out,keep_out,dkeep_out] = dmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS,id.KEEP,id.DKEEP);
     id.INFOG = inform;
     id.RINFOG = rinform;
     id.SOL = sol;
     id.INST = inst;
     id.SCHUR = schur;
     id.REDRHS = redrhs;
     id.PIVNUL_LIST = pivnul_list;
     id.SYM_PERM = sym_perm;
     id.UNS_PERM = uns_perm;
     id.TYPE = arithtype;
     id.ICNTL = icntl;
     id.CNTL = cntl;
     id.COLSCA = colsca_out;
     id.ROWSCA = rowsca_out;
     id.KEEP = keep_out;
     id.DKEEP = dkeep_out;
     return;
end

if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
end

if(id.TYPE ~= arithtype)
   disp('You are trying to call z/d version on a d/z instance');
   return;
end

[inform,rinform,sol,inst,schur,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl,colsca_out,rowsca_out,keep_out,dkeep_out] = dmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS,id.KEEP,id.DKEEP,mat);
id.INFOG = inform;
id.RINFOG = rinform;
id.SOL = sol;
id.INST = inst;
if(id.JOB == 2 | id.JOB == 4 | id.JOB == 6)
  if(id.SYM == 0)
        id.SCHUR = schur';
  else
        id.SCHUR = triu(schur)+tril(schur',-1);
  end
end
id.REDRHS = redrhs;
id.PIVNUL_LIST = pivnul_list;
id.SYM_PERM(sym_perm) = [1:size(mat,1)];
id.UNS_PERM = uns_perm;
id.ICNTL=icntl;
id.CNTL=cntl;
id.COLSCA=colsca_out;
id.ROWSCA=rowsca_out;
id.KEEP=keep_out;
id.DKEEP=dkeep_out;
