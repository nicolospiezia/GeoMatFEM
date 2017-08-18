function [AmI,ApI] = prodton(v,dof)
% -------------------------------------------------------------------------
% Compute tensor multiplication:
%
% {A (-) 1}_ijkl = A_il*I_jk = AmI
%
% {A (+) 1}_ijkl = A_jl*I_ik = ApI
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
% -------------------------------------------------------------------------

I = eye(3);
AI = zeros(3,3,3,3);
 
if dof == 2
    
    A   = [v(1)  v(3)   0
           v(3)  v(2)   0
           0     0      v(4)];
else
    
    A   = [v(1)  v(4)   v(6)
           v(4)  v(2)   v(5)
           v(6)  v(5)   v(3)];
end

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                
                AI(i,j,k,l) = A(i,l)*I(j,k);
               
            end
        end
    end
end 

if dof == 2
    % Allocate for a 2D analisys (non-symmetric)
    AmI =   [ AI(1,1,1,1)	AI(1,1,2,1)  AI(1,1,1,2) AI(1,1,2,2)
              AI(2,1,1,1)	AI(2,1,2,1)  AI(2,1,1,2) AI(2,1,2,2)
              AI(1,2,1,1)	AI(1,2,2,1)  AI(1,2,1,2) AI(1,2,2,2)
              AI(2,2,1,1)	AI(2,2,2,1)  AI(2,2,1,2) AI(2,2,2,2)];
else
    % Allocate for a 3D analisys 
    AmI=[AI(:,1,1,1) AI(:,1,2,1) AI(:,1,3,1) AI(:,1,1,2) AI(:,1,2,2) AI(:,1,3,2) AI(:,1,1,3) AI(:,1,2,3) AI(:,1,3,3)
         AI(:,2,1,1) AI(:,2,2,1) AI(:,2,3,1) AI(:,2,1,2) AI(:,2,2,2) AI(:,2,3,2) AI(:,2,1,3) AI(:,2,2,3) AI(:,2,3,3)
         AI(:,3,1,1) AI(:,3,2,1) AI(:,3,3,1) AI(:,3,1,2) AI(:,3,2,2) AI(:,3,3,2) AI(:,3,1,3) AI(:,3,2,3) AI(:,3,3,3)];
end
% -------------------------------------------------------------------------
if nargout == 2
    AI = zeros(3,3,3,3);

    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                 
                    AI(i,j,k,l) = A(j,l)*I(i,k);
               
                end
            end
        end
    end
    if dof == 2
        % Allocate for a 2D analisys (non-symmetric)
        ApI =   [ AI(1,1,1,1)	AI(1,1,2,1)  AI(1,1,1,2) AI(1,1,2,2)
                  AI(2,1,1,1)	AI(2,1,2,1)  AI(2,1,1,2) AI(2,1,2,2)
                  AI(1,2,1,1)	AI(1,2,2,1)  AI(1,2,1,2) AI(1,2,2,2)
                  AI(2,2,1,1)	AI(2,2,2,1)  AI(2,2,1,2) AI(2,2,2,2)];
    else
        % Allocate for a 3D analisys
ApI=[AI(:,1,1,1) AI(:,1,2,1) AI(:,1,3,1) AI(:,1,1,2) AI(:,1,2,2) AI(:,1,3,2) AI(:,1,1,3) AI(:,1,2,3) AI(:,1,3,3)
     AI(:,2,1,1) AI(:,2,2,1) AI(:,2,3,1) AI(:,2,1,2) AI(:,2,2,2) AI(:,2,3,2) AI(:,2,1,3) AI(:,2,2,3) AI(:,2,3,3)
     AI(:,3,1,1) AI(:,3,2,1) AI(:,3,3,1) AI(:,3,1,2) AI(:,3,2,2) AI(:,3,3,2) AI(:,3,1,3) AI(:,3,2,3) AI(:,3,3,3)];
    end
end
end

