function Dep = DEP_CAP21( epsev,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR)
%--------------------------------------------------------
% DEPtens: compute matrix Dep according to equatione 3.50 (Borja_Part III)
%
% Date: 29/10/2013
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------
% Inizialize matrices
[De,Dp] = deal(zeros(2));

% Compute parameter
Kp = pp*Pi0*((epstar/(epstar-epsevTR+epsev-epspvn))^(pp-1))*...
                (-epstar/(epstar-epsevTR+epsev-epspvn)^2);
            
KpTR = pp*Pi0*((epstar/(epstar-epsevTR+epsev-epspvn))^(pp-1))*...
    (epstar/(epstar-epsevTR+epsev-epspvn)^2);
            

% Compute derivative of Yield Function
dPF = 2*(B2)*(P-Pi);
dQF = 2*(A^2)*Q;
dPiF = (2*m^2*Pi+2*m*c0)*((P-Pi)^2-A^2)+B2*2*(P-Pi)*(-1);

% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;

% Compute matrix D

De = DE_EL(K,G);


% Compute parameters

d1 = De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF;
d2 = De(1,2)*dPF+De(2,2)*dQF;

e = d1*dPG+d2*dQG;

a1 = (d1+KpTR*dPiF)/e;
a2 = d2/e;

% Compute matrix Dp

Dp(1,1)= 1-a1*dPG;
Dp(1,2)= -a2*dPG;
Dp(2,1)= -a1*dQG;
Dp(2,2)= 1-a2*dQG;

% Compute matrix Dep

Dep = De*Dp;

end




