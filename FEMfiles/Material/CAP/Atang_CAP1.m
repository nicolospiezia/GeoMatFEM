function [Atang] = Atang_CAP1(K,G,m,mbar)

%ATANG Compute tangent for NR iteration for Cam Clay Return Mapping

% Inizialize matrices
[De] = deal(zeros(2));

Kp = 0;

% Compute derivative of Yield Function
dPF = -m;
dQF = 1;
dPiF = 0;

% Compute derivative of Plastic potential
dPG = -mbar;
dQG = 1;

% Compute matrix D

De = DE_EL(K,G);

% Compute matrix A

Atang = [                 1               0                              dPG
                          0               1                              dQG
          De(1,1)*dPF+De(2,1)*dQF+Kp*dPiF    De(1,2)*dPF+De(2,2)*dQF     0   ];
  
end




