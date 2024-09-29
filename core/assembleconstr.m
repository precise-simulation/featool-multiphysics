%ASSEMBLECONSTR Assemble integral constraints.
%
%   [ C, D ] = ASSEMBLECONSTR( PROB, ICUB ) Assemble integral constraints.
%
%   Assemble subdomain/boundary integral constraints specified in the
%   PROB.CONSTR field. Each constr entry have TYPE (intsubd/intbdr),
%   DVAR, INDEX (indices to subdomains/boundaries), and EXPR
%   (resulting integrand expression, must evaluate to a numeric
%   scalar) entries, for example
%
%       const(1).type = 'intsubd';   % Subdomain integral constraint.
%       const(1).dvar = 'u';         % Applied to dependent variable u.
%       const(1).index = [1,2];      % Applied to domains 1 and 2.
%       const(1).expr  = 2.3;        % Resulting integrand.

% Copyright 2013-2024 Precise Simulation, Ltd.
