function [tenspr,P,Q,epse,epspv,epsps,csi,aep,flag,fac,rmap] = tens(Ge,defepr,epspvn,epspsn,csin)
% -------------------------------------------------------------------------
% tens:
% Driver to compute the Kirchhoff principal tension for a general material
%
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% epspvn  : Plastic volumentric strain at step n
% epspsn  : Plastic deviatoric strain at step n
% csin    : Hardening parameter at step n
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]' 
%  epse  : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
% epspv  : Plastic volumentric strain at step n+1
% epsps  : Plastic deviatoric strain at step n+1
%  csi   : Hardening parameter at step n+1
%  aep   : Algorithmic stress-strain matrix in pr direction (a_AB).
%  flag  : Flag for plasticity.
%
% Date: 25/05/2015
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
% -------------------------------------------------------------------------

% Initialize vector 
tenspr = zeros(3,1);
epse = zeros(3,1);


% Compute volumetric trial strain
epsevTR = defepr(1,1)+defepr(2,1)+defepr(3,1);

% Compute deviatoric trial strain
epsedev(:,1) = defepr(:,1)-epsevTR/3;

epsesTR = sqrt(2/3)*norm(epsedev);

% Select material type
flagm=Ge(1);

switch flagm
      
    case 1 % Elastic
        epsev = epsevTR; 
        epses = epsesTR;
        [P,Q,Dep,flag,fac] = tens_EL(Ge,epsev,epses);
           
        csi = 0;
        epspv = 0;
        epsps = 0;
        
    case 2  % Drucker Prager
        cn = csin;    
        [P,Q,Dep,epsev,epses,epspv,epsps,c,flag,fac,rmap] = tens_DP(Ge,epsevTR,epsesTR,epspvn,epspsn,cn);
        csi = c;
    
    case 4  % Hyper Cam-Clay
        Pcn = csin;
        [P,Q,Dep,epsev,epses,epspv,epsps,Pc,flag,fac,rmap] = tens_CC(Ge,epsevTR,epsesTR,epspvn,epspsn,Pcn);
        csi = Pc;
        
    case 5  % Elastic Cap model
        Pin = csin;
        [P,Q,Dep,epsev,epses,epspv,epsps,Pi,flag,fac,rmap] = tens_CAP(Ge,epsevTR,epsesTR,epspvn,epspsn,Pin);
        csi = Pi;

        
    case 6  % Hyper Cap model
        Pin = csin;
        [P,Q,Dep,epsev,epses,epspv,epsps,Pi,flag,fac,rmap] = tens_CCDP(Ge,epsevTR,epsesTR,epspvn,epspsn,Pin);
        csi = Pi;
end

% Compute vector n

if norm(epsedev) == 0
n(:,1) = epsedev(:,1); 
else
n(:,1) = epsedev(:,1)/ norm(epsedev);  
end

% Compute principal Kirchhoff tension

tenspr(:,1) = P*ones(3,1)+ sqrt(2/3)*Q*n(:,1);

% Compute principal elastic strain

epse(:,1) = (1/3)*epsev*ones(3,1)+ sqrt(3/2)*epses*n(:,1);

% Compute algorithmic stress-strain in principal direction


if epsesTR == 0
    
%aep = [ K + 4*G/3  K - 2*G/3   K - 2*G/3 
%        K - 2*G/3  K + 4*G/3   K - 2*G/3
%        K - 2*G/3  K - 2*G/3   K + 4*G/3 ]; 
    
aep = (Dep(1,1)-2*Dep(2,2)/9)*ones(3,1)*ones(3,1)'+...
       sqrt(2/3)*Dep(1,2)*ones(3,1)*n(:,1)'+...
       sqrt(2/3)*Dep(2,1)*n(:,1)*ones(3,1)'+...
       ((2*Dep(2,2))/3)*eye(3); 

else
    
aep = (Dep(1,1)-2*Q/(9*epsesTR))*ones(3,1)*ones(3,1)'+...
       sqrt(2/3)*Dep(1,2)*ones(3,1)*n(:,1)'+...
       sqrt(2/3)*Dep(2,1)*n(:,1)*ones(3,1)'+...
      ((2*Q)/(3*epsesTR))*(eye(3)-n(:,1)*n(:,1)')+...
      (2/3)*Dep(2,2)*n(:,1)*n(:,1)';

end




% ------------------------------------------------------------ %
% --------------------- PLOT GRAFICO ------------------------- %
% ------------------------------------------------------------ %
%if ip == 1; % Plot solo del PG 1.

% P
% Q
    
    
% hold on
% 
% plot(P,Q,'bo')
% plot(Ptr,Qtr,'r*')
% axis([-700 0 0 600])
% axis equal
% xlabel('P')
% ylabel('Q')
% 
% % Draw ellipse initial
% xc = Pi;
% yc = 0;
% a = Ge(9);
% b = sqrt(m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2); % 
% angle = 0;
% ell = calculateEllipse(xc, yc, a, b, angle,200);
% plot(ell(:,1), ell(:,2), '-'), axis([-700 0 0 600])
% 
% % Draw ellipse current
% xc = Pi0;
% yc = 0;
% a = Ge(9);
% b = sqrt(m^2*Pi0^2-m^2*A^2+2*m*c0*Pi0+c0^2); % 
% angle = 0;
% ell = calculateEllipse(xc, yc, a, b, angle,200);
% plot(ell(:,1), ell(:,2), ':b'), axis([-700 0 0 600])
% 
% 
% % Draw tangent line
% pstar = Pi-(m*a*b)/(sqrt((b^4/a^2)+m^2*b^2));
% qstar = Ge(5)*pstar+Ge(10);
% pl1 = [pstar pstar];
% ql2 = [0 qstar];
% plot(pl1, ql2, ':g' )
% 
% % Draw line
% pl = 0:-1:pstar;
% ql = Ge(5)*pl+Ge(10);
% plot(pl, ql, '--r' ), axis([-700 0 0 600])
% %   
%end
% % % ------------------------------------------------------------ %

end  

