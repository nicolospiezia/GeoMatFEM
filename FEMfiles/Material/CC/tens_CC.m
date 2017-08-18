function [P,Q,Dep,epsev,epses,epspv,epsps,Pc,flag,fac,NORMErec] = tens_CC(Ge,epsevTR,epsesTR,epspvn,epspsn,Pcn)
%--------------------------------------------------------
% tensVM: 
%   Compute the Kirchhoff principal tension according to CC yield criterion.
%
% Syntax:
%   [tenspr,epse,Pc,aep,flag] = tensCC(Ge,defepr,Pcn)
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% Pcn     : Preconsolidation pressure at step n.
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]' 
%   epse : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
%   Pc   : Preconsolidation pressure at step n+1.
%  aep   : Algorithmic stress-strain matrix in pr direction (a_AB).
%  flag  : Flag for plasticity.
%
% Reference:
% Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
% infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
%
% Date:
%   Version 1.0   25.10.13
%
% Created by: Nicolò Spiezia
%--------------------------------------------------------
%load gong.mat;
imax = 15;
toll = 10e-7; % 6 9
NORMErec = zeros(imax,1);

% Set isotropic elasto-plastic parameter
mu0  = Ge(3);  
alfa = Ge(4);
kappa = Ge(5);
lambda = Ge(6);
THETA = 1/(lambda-kappa);
M = Ge(7);
P0 = Ge(8);
epsev0 = Ge(10);

% Compute trial stress invariants
[Ptr,Qtr] = PQ_HY( epsevTR, epsesTR, P0, alfa, kappa, epsev0, mu0);

% Check for plasticity    
Ftr = (Qtr/M)^2+Ptr*(Ptr-Pcn);
    
if Ftr > toll % PLASTIC STEP    
   flag = 1;
   
        % Solve NR system
        
            % Inizialize variables
            x = zeros(3,1);
            x(1,1) = epsevTR;   % epsev
            x(2,1) = epsesTR;   % epses
            x(3,1) = 0;         % dgamma
                       
            Pc = Pcn;
            P = Ptr;
            Q = Qtr;
            
            for iter = 1:imax
        
                % evaluate residual            
                r = [x(1,1) - epsevTR + x(3,1)*(2*P-Pc);
                     x(2,1) - epsesTR + x(3,1)*(2*Q/M^2);
                     (Q/M)^2 + P*(P-Pc)];

                 
                if iter == 1
                r0 = norm(r);
                end
                      
                norme = norm(r)/r0;
                NORMErec(iter,1)=norme; 
                
                % check for convergence
                if norme < toll    
                    break
                else

                    % evaluate tangent matrix for NR iteration
                    A = Atang_CC(x(1,1),x(2,1),x(3,1),P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M);

                    % solve for displacement increment
                    dx = - (A\r);

                    x = x + dx;                 
                    
                    % Update 
                    
                    [P,Q] = PQ_HY(x(1,1),x(2,1),P0,alfa,kappa,epsev0,mu0);
                    
                    epspv = epsevTR-x(1)+epspvn;
                    epsps = epsesTR-x(2)+epspsn;
                    
                    Pc = Pcn*exp(-THETA*(epsevTR-x(1,1)));
                    
                        if iter == imax
                        fprintf('\n No convergence RM \n');
                        %sound(y);
                        break
                        end
                
                end
            end
            
        %Update variable
        epsev = x(1,1);
        epses = x(2,1);
        dgamma = x(3,1);
        
        % Compute matrix Dep (Eq. 3.50)
        
        Dep = DEP_CC(epsev,epses,dgamma,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M);
        
       
else % ELASTIC STEP
     flag = 0;
        
        %Update variable
        epsev = epsevTR;
        epses = epsesTR;
        
        epspv = epspvn;
        epsps = epspsn;
        
        Pc = Pcn;
        P = Ptr;
        Q = Qtr;
        
        % Compute matrix De 
        
        De = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0);
        
        Dep = De;
        
end

OMEGA = -(epsev-epsev0)/kappa;
fac = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));  % Factor for computing alfa tensor in case of NaN


  


end

