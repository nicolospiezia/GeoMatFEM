function [A] = Atang_CC( epsev,epses,dgamma,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M )
%--------------------------------------------------------
% ATANG_CC Compute tangent for NR iteration for Cam Clay Return Mapping
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
A = zeros(3);
[De, H, G] = deal(zeros(2));

% Compute parameter
OMEGA = -(epsev-epsev0)/kappa;
mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));
THETA = 1/(lambda-kappa);


% Compute matrix D

De = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0);

% Compute matrix H

H(1,1) = 2;
H(2,2) = 2/(M^2);
H(1,2) = 0;
H(2,1) = 0;


% Compute matrix G

G = H*De;

% Compute matrix A

A(1,1) = 1+dgamma*(G(1,1)-THETA*Pc);
A(1,2) = dgamma*G(1,2);
A(1,3) = 2*P-Pc;

A(2,1) = dgamma*G(2,1);
A(2,2) = 1+dgamma*G(2,2);
A(2,3) = 2*Q/(M^2);

A(3,1) = De(1,1)*(2*P-Pc)+De(2,1)*(2*Q/(M^2))+THETA*Pc*(-P);
A(3,2) = De(1,2)*(2*P-Pc)+De(2,2)*(2*Q/(M^2));
A(3,3) = 0;

end




