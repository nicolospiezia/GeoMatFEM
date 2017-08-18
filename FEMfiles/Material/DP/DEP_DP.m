function Dep = DEP_DP(K,G,m,mbar,H)
%--------------------------------------------------------
% DEPtens: compute matrix Dep for Drucker-Prager 
%
% Date: 14/01/2014
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------
% Inizialize matrices
[De,Dp] = deal(zeros(2));


% Compute derivative of Yield Function
dPF = -m;
dQF = 1;
dcF = -1;
Kp = -H;
KpTr = H;

% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;

% Compute matrix D

De = DE_EL(K,G);

% Compute parameters

d1 = De(1,1)*dPF+De(2,1)*dQF+Kp*dcF;
d2 = De(1,2)*dPF+De(2,2)*dQF;

e = d1*dPG+d2*dQG;

a1 = (d1+KpTr*dcF)/e;
a2 = d2/e;

% Compute matrix Dp

Dp(1,1) = 1-a1*dPG;
Dp(1,2) = -a2*dPG;
Dp(2,1) = -a1*dQG;
Dp(2,2) = 1-a2*dQG;

% Compute matrix Dep

Dep = De*Dp;

end




