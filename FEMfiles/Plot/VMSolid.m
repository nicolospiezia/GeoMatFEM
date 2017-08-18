function f = VMSolid(S)
%--------------------------------------------------------  
% VMHex:
%   Calculates the value of the von Mises
%   yield function(equivalent stress).
%
% Syntax:
%   f = VMHex(S)
%
% Input:
%   S   : Stress array.  
%
% Output:
%   f   : Value of yield function.
%
% Date:
%   Version 3.0   31.10.14
%--------------------------------------------------------

% System parameters
nel = size(S,3);
nip = size(S,2);

% Initialize f
f = zeros(nel,nip);

for i = 1:nel

    for j = 1:nip
        
        % Extract stresses
        Ss = S(:,j,i);
        
        if size(S,1)==4
            
            % Von Mises (equivalent stress)        
            %f(i,j) = sqrt(Ss(1)^2+Ss(2)^2-Ss(1)*Ss(2)+3*Ss(3)^2); spie
            f(i,j) = sqrt(Ss(1)^2+Ss(2)^2+Ss(4)^2-Ss(1)*Ss(2)-Ss(1)*Ss(4)-Ss(2)*Ss(4)+3*Ss(3)^2);
        else
        
            % Intermidiate variables
            S12 = Ss(1) - Ss(2);
            S13 = Ss(1) - Ss(3);
            S23 = Ss(2) - Ss(3);
            t23 = Ss(4);
            t31 = Ss(5);
            t12 = Ss(6);
               
            % Octehedral shear stress
            t0 = sqrt( S12^2 + S13^2 + S23^2 + ...
                     6*t12^2 + 6*t31^2 + 6*t23^2 ); 
        
            % Von Mises (equivalent stress)        
            f(i,j) = 1/sqrt(2)*t0;
        end
    end
end
