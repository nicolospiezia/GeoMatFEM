function [ brec,PLrec ] = InizbPL(nip,dof,ne,nelem,nonlin,T,G)
%INIZBPL Summary of this function goes here
%   Detailed explanation goes here

brec=zeros(2*dof,nip^dof,nelem);
PLrec=zeros(1,nip^dof,nelem);

if dof == 2
    m = [1 1 0 1]';
else
    m = [1 1 1 0]';
end

for i=1:nelem;
    for j=1:nip^dof 
        
        if nonlin == 1
        % Initial elastic left Cauchy-Green tensor
        brec(1:4,j,i) = m;
        end
        
        % Plastic variable record
        flag = G(T(i,ne+1),1);
    
        switch flag
            
            case 2 % Drucker-Prager
                PLrec (1,j,i) = G(T(i,ne+1),6); % Allocate c0
            
            case 4 % Hyper Cam-Clay
                PLrec (1,j,i) = G(T(i,ne+1),9); % Allocate Pc0
                
            case 5 % Elastic Cap model
                PLrec (1,j,i) = G(T(i,ne+1),6); % Allocate Pi
                
            case 6 % Hyper Cap Model
                PLrec (1,j,i,1) = G(T(i,ne+1),8); % Allocate Pi
        
        end    
    end
end
end

