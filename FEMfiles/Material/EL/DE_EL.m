function DE = DE_EL(K,G)
% DEtens Compute the elastic tensor De (3.33a Part III) 
% for isotropic linear elastic material


% Inizialize matrices
De = zeros(2);


% Compute elastic tensor De (3.33a Part III)
DE(1,1) = K;             
DE(2,2)= 3*G;

end




