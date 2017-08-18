function [alfamat,alfatens] = alfa(aep,tenspr,lambar,avett,dof,nonlin,fac)
%--------------------------------------------------------
% Compute tangent operator and store in vector form.
% 
% Date: 31/10/2014
%   Version 3.0    
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Compute the first term:

alfa1 = zeros(3,3,3,3);
alfa2 = zeros(3,3,3,3);
%alfavet=zeros(81,1);
    
for A=1:3
    for B=1:3   
        
        % Compute the first term:
        
        alfa1 = alfa1+aep(A,B)*dyadic(avett(:,A)*avett(:,A)',avett(:,B)*avett(:,B)');
        
        % Compute the second term:
    
        if B ~= A
            
            if lambar(B)~=lambar(A)
            
                fac =  ((tenspr(B)-tenspr(A))/(lambar(B)-lambar(A)));
                
            end
            %lambar(B)==lambar(A)
            %fac = fac = G        -> ELASTIC  
            %fac = fac = mu0+...  -> HYPER-ELASTIC
            
            if nonlin==0
                alfa2 = alfa2+fac*...
                (dyadic(avett(:,A)*avett(:,B)',avett(:,A)*avett(:,B)')+...
                 dyadic(avett(:,A)*avett(:,B)',avett(:,B)*avett(:,A)'))/2;
            else
                alfa2 = alfa2+fac*...
                (lambar(B)*dyadic(avett(:,A)*avett(:,B)',avett(:,A)*avett(:,B)')+...
                 lambar(A)*dyadic(avett(:,A)*avett(:,B)',avett(:,B)*avett(:,A)'));
                
            end
        end
    
    end
end
%
% Total alfa tensor 
alfatens = alfa1+alfa2;

if dof==2
    
% Allocate alfa as vector for 2D plain strain (non-symmetric)

alfamat=[alfatens(1,1,1,1)	alfatens(1,1,2,1)  alfatens(1,1,1,2) alfatens(1,1,2,2)
         alfatens(2,1,1,1)	alfatens(2,1,2,1)  alfatens(2,1,1,2) alfatens(2,1,2,2)
         alfatens(1,2,1,1)	alfatens(1,2,2,1)  alfatens(1,2,1,2) alfatens(1,2,2,2)
         alfatens(2,2,1,1)	alfatens(2,2,2,1)  alfatens(2,2,1,2) alfatens(2,2,2,2)];

%alfavet = reshape(alfamat,16,1);
else

% Allocate alfa as vectos for 3D solid hexaedral element 

%alfamat=[alfatens(1,1,1,1) alfatens(1,1,2,1) alfatens(1,1,3,1) alfatens(1,1,1,2) alfatens(1,1,2,2) alfatens(1,1,3,2) alfatens(1,1,1,3) alfatens(1,1,2,3) alfatens(1,1,3,3)
%         alfatens(2,1,1,1) alfatens(2,1,2,1) alfatens(2,1,3,1) alfatens(2,1,1,2) alfatens(2,1,2,2) alfatens(2,1,3,2) alfatens(2,1,1,3) alfatens(2,1,2,3) alfatens(2,1,3,3)
%         alfatens(3,1,1,1) alfatens(3,1,2,1) alfatens(3,1,3,1) alfatens(3,1,1,2) alfatens(3,1,2,2) alfatens(3,1,3,2) alfatens(3,1,1,3) alfatens(3,1,2,3) alfatens(3,1,3,3)
%         alfatens(1,2,1,1) alfatens(1,2,2,1) alfatens(1,2,3,1) alfatens(1,2,1,2) alfatens(1,2,2,2) alfatens(1,2,3,2) alfatens(1,2,1,3) alfatens(1,2,2,3) alfatens(1,2,3,3)
%         alfatens(2,2,1,1) alfatens(2,2,2,1) alfatens(2,2,3,1) alfatens(2,2,1,2) alfatens(2,2,2,2) alfatens(2,2,3,2) alfatens(2,2,1,3) alfatens(2,2,2,3) alfatens(2,2,3,3)
%         alfatens(3,2,1,1) alfatens(3,2,2,1) alfatens(3,2,3,1) alfatens(3,2,1,2) alfatens(3,2,2,2) alfatens(3,2,3,2) alfatens(3,2,1,3) alfatens(3,2,2,3) alfatens(3,2,3,3)
%         alfatens(1,3,1,1) alfatens(1,3,2,1) alfatens(1,3,3,1) alfatens(1,3,1,2) alfatens(1,3,2,2) alfatens(1,3,3,2) alfatens(1,3,1,3) alfatens(1,3,2,3) alfatens(1,3,3,3)
%         alfatens(2,3,1,1) alfatens(2,3,2,1) alfatens(2,3,3,1) alfatens(2,3,1,2) alfatens(2,3,2,2) alfatens(2,3,3,2) alfatens(2,3,1,3) alfatens(2,3,2,3) alfatens(2,3,3,3)
%         alfatens(3,3,1,1) alfatens(3,3,2,1) alfatens(3,3,3,1) alfatens(3,3,1,2) alfatens(3,3,2,2) alfatens(3,3,3,2) alfatens(3,3,1,3) alfatens(3,3,2,3) alfatens(3,3,3,3)];

alfamat=[alfatens(:,1,1,1) alfatens(:,1,2,1) alfatens(:,1,3,1) alfatens(:,1,1,2) alfatens(:,1,2,2) alfatens(:,1,3,2) alfatens(:,1,1,3) alfatens(:,1,2,3) alfatens(:,1,3,3)
         alfatens(:,2,1,1) alfatens(:,2,2,1) alfatens(:,2,3,1) alfatens(:,2,1,2) alfatens(:,2,2,2) alfatens(:,2,3,2) alfatens(:,2,1,3) alfatens(:,2,2,3) alfatens(:,2,3,3)
         alfatens(:,3,1,1) alfatens(:,3,2,1) alfatens(:,3,3,1) alfatens(:,3,1,2) alfatens(:,3,2,2) alfatens(:,3,3,2) alfatens(:,3,1,3) alfatens(:,3,2,3) alfatens(:,3,3,3)];

%alfavet = reshape(alfamat,81,1);
%--------------------------------------------------------------------------
%m=0;
%for i=1:3
%    for j=1:3
%        for k=1:3
%            for l=1:3
%                m=m+1;
%                alfavet(m)=alfatens(l,k,j,i);
%            end
%        end
%    end
%end
end
end

               