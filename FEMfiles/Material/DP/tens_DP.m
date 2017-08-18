function [P,Q,Dep,epsev,epses,epspv,epsps,c,flag,fac,NORMErec] = tens_DP(Ge,epsevTR,epsesTR,epspvn,epspsn,cn)
%--------------------------------------------------------
% tensVM: 
%   Compute the Kirchhoff principal tension according to 
%   associative or non associative Drucker-Prager yield criterion.
%
% Syntax:
%   [tenspr,P,Q,epse,epspv,aep,flag] = tensDP(Ge,defepr,epspvn)
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% epspvn  : Plastic volumetric deformation at step n.
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]'
% P      : Mean normal Kirchhoff stress
% Q      : Deviatoric Kirchhoff stress      
% epse   : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
% epspv  : Plastic volumetric deformation 
% aep    : Algorithmic stress-strain matrix in principal direction (a_AB)
% flag   : Flag for plasticity
%
% Date:
%   Version 1.0   25.09.14
%--------------------------------------------------------

% Iteration parameters
imax = 15;
toll = 10e-7;

% Initialize vector 
NORMErec = zeros(imax,1);

% Set isotropic elasto-plastic parameter

%  3	 4	  5	    6	  7	  8
%  E    ni   m    c0    mbar  H

E = Ge(3);  
ni = Ge(4);
m = Ge(5);
c0 = Ge(6);
mbar = Ge(7);
H = Ge(8); % Hardening modulus

G = E/(2*(1+ni));
K = E/(3*(1-2*ni));

 
% Compute TRial stress    

[Ptr,Qtr] = PQ_EL(epsevTR, epsesTR,K,G);

flag = 0;

% DP trial yield surface
Ftr = Qtr-m*Ptr-cn;
        
          
 % Check for plasticity

if Ftr > toll 
flag = 1;    
        
        % !------------ Drucker Prager RM --------------!
            % Inizialize variables
            x = zeros(3,1);
            x(1,1) = epsevTR;   % epsev
            x(2,1) = epsesTR;   % epses
            x(3,1) = 0;         % dgamma
            
            c = cn;
            P = Ptr;
            Q = Qtr;
            
            for iter = 1:imax
      
                % evaluate residual                
                
                F = Q-m*P-c;

                dPG = -mbar;
                dQG = 1;

                r = [x(1) - epsevTR+x(3)*dPG
                     x(2) - epsesTR+x(3)*dQG
                           F                ];
                
                   if iter == 1
                   r0 = norm(r);
                   end
                       
                norme = norm(r)/r0;
                NORMErec(iter,1)=norme; 

                % check for convergence
                if norme < toll
                    break
                else

                    % evaluate tangent matrix
                    Atang = Atang_DP(K,G,m,mbar,H);
                
                    % solve for displacement increment
                    dx = - (Atang\r);

                    x = x + dx;                 
                    
                    % Update 
                    
                    [P,Q] = PQ_EL(x(1,1), x(2,1),K,G);
                    
                    epspv = epsevTR-x(1)+epspvn;
                    epsps = epsesTR-x(2)+epspsn;
                    
                    c = c0+H*epspv;
 
                end
            end
            
                if iter == imax
                fprintf('No convergence of RM')
                end
           
       
        %Update variable
        epsev = x(1,1);
        epses = x(2,1);
        dgamma = x(3,1);
        
        % Compute matrix Dep for DP non-associative
        
        Dep = DEP_DP(K,G,m,mbar,H);
        % !---------------------------------------------!
        
                
       % Check for apex RM
        if Q < 0
        fprintf('Apex region RM')
        end
      
else   % ELASTIC STEP
     
    %Update variable
    epsev = epsevTR;
    epses = epsesTR;

    P = Ptr;
    Q = Qtr;
    epspv = epspvn;
    epsps = epspsn;
    c = cn;

    % Compute matrix De 

    De = DE_EL(K,G);

    Dep = De;
        
end
   
   fac = G;


end

