function Dep = DEP_CC( epsev,epses,dgamma,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M )
%--------------------------------------------------------
% DEP_CC: compute matrix Dep according to equatione 3.50 
%
% Reference:
% Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
% infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
%
% Date: 25/10/2013
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------

% Inizialize matrices
[H,b,Dp] = deal(zeros(2));

% Compute parameter
OMEGA = -(epsev-epsev0)/kappa;
mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));
THETA = 1/(lambda-kappa);

% Compute matrix De

De = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0);

% Compute matrix H

H(1,1) = 2;
H(2,2) = 2/(M^2);
H(1,2) = 0;
H(2,1) = 0;

% Compute matrix G

G = H*De;

% Compute matrix b

b(1,1) = 1+dgamma*(G(1,1)-THETA*Pc);
b(1,2) = dgamma*G(1,2);
b(2,1) = dgamma*G(2,1);
b(2,2) = 1+dgamma*G(2,2);

% Compute parameters

c1 = 1-dgamma*(-THETA*Pc)*(-1); 
c2 = -dgamma*(-THETA*Pc)*0;

d1 = De(1,1)*(2*P-Pc)+De(2,1)*(2*Q/(M^2))+THETA*Pc*(-P);
d2 = De(1,2)*(2*P-Pc)+De(2,2)*(2*Q/(M^2));

e = d1*(b(2,2)*(2*P-Pc)-b(1,2)*(2*Q/(M^2)))+d2*(b(1,1)*(2*Q/(M^2))-b(2,1)*(2*P-Pc));

a1 = (d1*(b(2,2)*c1-b(1,2)*c2)+d2*(b(1,1)*c2-b(2,1)*c1)+det(b)*(-THETA*Pc)*(-P))/e;
a2 = sqrt(2/3)*(d2*b(1,1)-d1*b(1,2))/e;

% Compute matrix Dp

Dp(1,1)= b(2,2)*(c1-a1*(2*P-Pc))-b(1,2)*(c2-a1*(2*Q/(M^2)));
Dp(1,2)= b(1,2)*(-1+sqrt(3/2)*a2*(2*Q/(M^2)))-sqrt(3/2)*b(2,2)*a2*(2*P-Pc);
Dp(2,1)= b(1,1)*(c2-a1*(2*Q/(M^2)))-b(2,1)*(c1-a1*(2*P-Pc));
Dp(2,2)= b(1,1)*(1-sqrt(3/2)*a2*(2*Q/(M^2)))+sqrt(3/2)*b(2,1)*a2*(2*P-Pc);

Dp = Dp/det(b);

% Compute matrix Dep

Dep = De*Dp;

end




