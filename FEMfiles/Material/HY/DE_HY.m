function DE = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0)
%DE_HY Compute the elastic tensor De (3.33a Part III)


% Inizialize matrices
De = zeros(2);

% Compute parameter
OMEGA = -(epsev-epsev0)/kappa;
mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));

% Compute elastic tensor De (3.33a Part III)
DE(1,1) = -P/kappa;             
DE(2,2)= 3*mue;
DE(1,2)= (3*P0*alfa*epses/kappa)*exp(OMEGA);
DE(2,1)= (3*P0*alfa*epses/kappa)*exp(OMEGA);

end




