function [Atang] = Atang_CCDP2( epsev,epses,dgamma,P,Q,Pi,kappa,mu0,alfa,P0,epsev0,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR)

%ATANG Compute tangent for NR iteration for Cam Clay Return Mapping

% Inizialize matrices
[De] = deal(zeros(2));

% Compute parameter
OMEGA = -(epsev-epsev0)/kappa;
mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));

Kp = pp*Pi0*((epstar/(epstar-epsevTR+epsev-epspvn))^(pp-1))*...
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

De(1,1) = -P/kappa;             
De(2,2)= 3*mue;
De(1,2)= (3*P0*alfa*epses/kappa)*exp(OMEGA);
De(2,1)= (3*P0*alfa*epses/kappa)*exp(OMEGA);

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




