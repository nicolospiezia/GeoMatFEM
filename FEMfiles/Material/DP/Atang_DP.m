function [Atang] = Atang_DP(K,G,m,mbar,H)

%ATANG Compute tangent for NR iteration for Drucker-Prager Return Mapping

% Inizialize matrices
[De] = deal(zeros(2));

Kp = -H;

% Compute derivative of Yield Function
dPF = -m;
dQF = 1;
dcF = -1;

% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;


% Compute matrix D

De = DE_EL(K,G);
 
% Compute matrix A

Atang = [                 1               0                              dPG
                          0               1                              dQG
          De(1,1)*dPF+De(2,1)*dQF+Kp*dcF    De(1,2)*dPF+De(2,2)*dQF     0   ];
  
end




