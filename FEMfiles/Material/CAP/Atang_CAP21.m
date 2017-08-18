function [Atang] = Atang_CAP21( epsev,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR)

%ATANG Compute tangent for NR iteration for Cam Clay Return Mapping

% Inizialize matrices
[De] = deal(zeros(2));

% Compute parameter

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

De = DE_EL(K,G);

% Compute matrix A

Atang = [   1                            0                              dPG
            0                            1                              dQG
          De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF    De(1,2)*dPF+De(2,2)*dQF     0   ];
  
end




