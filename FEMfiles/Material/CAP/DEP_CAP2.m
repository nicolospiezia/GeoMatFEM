function Dep = DEP_CAP2( epsev,dgamma,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR)
%--------------------------------------------------------
% DEPtens: compute matrix Dep for CC for rocks
%
% Date: 29/10/2013
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------
% Inizialize matrices
[De,b,Dp] = deal(zeros(2));

% Compute parameter
Kp = pp*Pi0*((epstar/(epstar-epsevTR+epsev-epspvn))^(pp-1))*...
                (-epstar/(epstar-epsevTR+epsev-epspvn)^2);
            
KpTR = pp*Pi0*((epstar/(epstar-epsevTR+epsev-epspvn))^(pp-1))*...
                (epstar/(epstar-epsevTR+epsev-epspvn)^2);

% Compute derivative of Yield Function
dPF = 2*(B2)*(P-Pi);
dQF = 2*(A^2)*Q;
dPiF = (2*m^2*Pi+2*m*c0)*((P-Pi)^2-A^2)+B2*2*(P-Pi)*(-1);

dPPF = 2*B2;
dQQF = 2*(A^2);
dPQF = 0;

dPPiF = 2*((2*m^2*Pi+2*m*c0)*(P-Pi)-B2);
dQPiF = 0;

% Compute matrix D

De = DE_EL(K,G);

% Compute matrix H

H = [ dPPF  dPQF
      dPQF  dQQF];


% Compute matrix G

G = H*De;

% Compute matrix b

b(1,1) = 1+dgamma*(G(1,1)+Kp*dPPiF);
b(1,2) = dgamma*G(1,2);
b(2,1) = dgamma*(G(2,1)+Kp*dQPiF);
b(2,2) = 1+dgamma*G(2,2);

% Compute parameters

c1 = 1-dgamma*KpTR*dPPiF; 
c2 = -dgamma*KpTR*dQPiF;

d1 = De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF;
d2 = De(1,2)*dPF+De(2,2)*dQF;

e = d1*(b(2,2)*dPF-b(1,2)*dQF)+d2*(b(1,1)*dQF-b(2,1)*dPF);

a1 = (d1*(b(2,2)*c1-b(1,2)*c2)+d2*(b(1,1)*c2-b(2,1)*c1)+det(b)*KpTR*dPiF)/e;
a2 = sqrt(2/3)*(d2*b(1,1)-d1*b(1,2))/e;

% Compute matrix Dp

Dp(1,1)= b(2,2)*(c1-a1*dPF)-b(1,2)*(c2-a1*dQF);
Dp(1,2)= b(1,2)*(-1+sqrt(3/2)*a2*dQF)-sqrt(3/2)*b(2,2)*a2*dPF;
Dp(2,1)= b(1,1)*(c2-a1*dQF)-b(2,1)*(c1-a1*dPF);
Dp(2,2)= b(1,1)*(1-sqrt(3/2)*a2*dQF)+sqrt(3/2)*b(2,1)*a2*dPF;

Dp = Dp/det(b);

% Compute matrix Dep

Dep = De*Dp;

end




