function Dep = DEP__CCDP1( epsev,epses,P,kappa,mu0,alfa,P0,epsev0,m,mbar)
%--------------------------------------------------------
% DEPtens: compute matrix Dep for DP non-associative
%
% Date: 14/01/2014
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------
% Inizialize matrices
[De,Dp] = deal(zeros(2));

% Compute parameter
OMEGA = -(epsev-epsev0)/kappa;
mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));

% Compute derivative of Yield Function
dPF = -m;
dQF = 1;

% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;

% Compute matrix D

De = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0);

% Compute parameters

d1 = De(1,1)*dPF+De(2,1)*dQF;
d2 = De(1,2)*dPF+De(2,2)*dQF;

e = d1*dPG+d2*dQG;

a1 = d1/e;
a2 = d2/e;

% Compute matrix Dp

Dp(1,1)= 1-a1*dPG;
Dp(1,2)= -a2*dPG;
Dp(2,1)= -a1*dQG;
Dp(2,2)= 1-a2*dQG;

% Compute matrix Dep

Dep = De*Dp;

end




