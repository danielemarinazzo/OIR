%% Adapted from SSGC Toolbox of Lionel Barnett
% L. Barnett and A. K. Seth, Granger causality for state-space models, Phys.Rev. E 91(4) Rapid Communication, 2015.

function [A,C,K,rho] = oir_ar2iss(Am)
    
% Return innovations form state space model parameters corresponding to a vector autoregressive model.
%
% Am=[A(1)...A(p)]: M*pM matrix of the MVAR model coefficients
%
% A,C,K - innovations form state space parameters
%
% rho   - AR spectral norm
%
%Note that rho >= 1 indicates an unstable AR process: rho > 1 is explosive, rho close to 1 may be unit-root.

M=size(Am,1);
p = size(Am,2)/M; % p is VAR model order

pM1 = (p-1)*M;

C = Am;
A = [C; eye(pM1) zeros(pM1,M)];
K = [eye(M); zeros(pM1,M)];

if nargout > 3
    rho = max(abs(eig(A,'nobalance')));
end
