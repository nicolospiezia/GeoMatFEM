function [Atang] = Atang_CAP2( epsev,dgamma,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR)

%ATANG Compute tangent for NR iteration for Cam Clay Return Mapping

% Inizialize matrices
[De] = deal(zeros(2));

% Compute parameter
epspv = epsevTR-epsev+epspvn;
Kp = pp*Pi0*((epstar/(epstar-epspv))^(pp-1))*...
                (-epstar/(epstar-epsevTR+epsev-epspvn)^2);

% Compute derivative of Yield Function
dPF = 2*B2*(P-Pi);
dQF = 2*(A^2)*Q;
dPiF = (2*Pi*m^2+2*m*c0)*((P-Pi)^2-A^2)+B2*2*(P-Pi)*(-1);

dPPF = 2*B2;
dQQF = 2*(A^2);
dPQF = 0;

dPPiF = 2*((2*Pi*m^2+2*m*c0)*(P-Pi)-B2);
dQPiF = 0;

% Compute matrix D

De = DE_EL(K,G);

% Compute matrix H

H = [ dPPF  dPQF
      dPQF  dQQF];
  
% Compute matrix G

G = H*De;

% Compute matrix A

Atang = [   1+dgamma*(G(1,1)+Kp*dPPiF)             dgamma*G(1,2)        dPF
             dgamma*(G(2,1)+Kp*dQPiF)             1+dgamma*G(2,2)       dQF
          De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF    De(1,2)*dPF+De(2,2)*dQF     0   ];
  
end




