function [Atang] = Atang_CCDP21( epsev,epses,P,Q,Pi,kappa,mu0,alfa,P0,epsev0,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR)

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


% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;

% Compute matrix D

De(1,1) = -P/kappa;             
De(2,2)= 3*mue;
De(1,2)= (3*P0*alfa*epses/kappa)*exp(OMEGA);
De(2,1)= (3*P0*alfa*epses/kappa)*exp(OMEGA);

% Compute matrix A

Atang = [   1                            0                              dPG
            0                            1                              dQG
          De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF    De(1,2)*dPF+De(2,2)*dQF     0   ];
  
end




