function [P,Q,Dep,epsev,epses,epspv,epsps,Pi,flag,fac,NORMErec] = tens_CCDP(Ge,epsevTR,epsesTR,epspvn,epspsn,Pin)

%--------------------------------------------------------
% tensVM: 
%   Compute the Kirchhoff principal tension according to CC yield criterion.
%
% Syntax:
%   [tenspr,PL,flag] = tensCC(Ge,defepr,Pcn)
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% epbar   : Preconsolidation pressure at step n.
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]' 
%   epse : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
%   Pc   : Preconsolidation pressure at step n+1
%  aep   : Algorithmic stress-strain matrix in pr direction (a_AB)
%  flag  : Flag for plasticity
%
% Date:
%   Version 1.0   25.10.13
%--------------------------------------------------------

% Initialize vector
imax = 20;
F_toll = 10e-4; % Attenzione: la tolleranza deve essere coerente 
                % con la grandezza di PQ del problema

tenspr = zeros(3,1);
epse = zeros(3,1);
NORMErec = zeros(imax,1);

% Set isotropic elasto-plastic parameter

% 2 	3	  4	     5	6	7	 8	     9	10	11	12  13
% rho  mu0  alfa  kappa  m  P0  Pi0  epsev0  A  c0  pp  epstar

mu0  = Ge(3);  
alfa = Ge(4);
kappa = Ge(5);
m = Ge(6);
P0 = Ge(7);
Pi0 = Ge(8);
epsev0 = Ge(9);
A = Ge(10);
c0 = Ge(11);
pp = Ge(12);
epstar = Ge(13);
mbar = Ge(14);


% Compute TRial stress    

[Ptr,Qtr] = PQ_HY( epsevTR, epsesTR, P0, alfa, kappa, epsev0, mu0);

flag = 0;

% Compute trial yield surface 

        % CC yield surface
        B2 = m^2*Pin^2-m^2*A^2+2*m*c0*Pin+c0^2;
        F1tr = B2*(Ptr-Pin)^2+(A^2)*(Qtr^2)-(A^2)*B2;

        % DP yield surface
        F2tr = Qtr-m*Ptr-c0;
        
        % Intersection CC/DP
        a = A;
        b = sqrt(m^2*Pin^2-m^2*A^2+2*m*c0*Pin+c0^2); 
        pstar = Pin-(m*a*b)/(sqrt((b^4/a^2)+m^2*b^2));
           
 % Check for plasticity

if F1tr > F_toll 
  
  if F2tr > F_toll 
     flag = 1;    
     surf = 2; % Drucker prager
        
    [P,Q,epsev,epses,epspv,epsps,Pi,Dep] = rm_DP(epsevTR,epsesTR,Ptr,Qtr,...
                                           epspvn,epspsn,kappa,mu0,alfa,P0,epsev0,m,c0,mbar,Pin);
                                       
        
                
               % Check for apex RM
                if Q < 0
                fprintf('Apex region RM')
                end
       
        
                % Check for correct surface
                B2 = m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2;
                F1 = B2*(P-Pi)^2+(A^2)*(Q^2)-(A^2)*B2;



                if F1 > F_toll && P < pstar
                   surf = 1; % Goto Cam Clay surface
                end
                
        else
            if Ptr<pstar
            flag = 1;
            surf = 1;   
            end         
      
        end

   end
       
  
  if flag ==1 && surf == 1
        
 [P,Q,epsev,epses,epspv,epsps,Pi,Dep,flag] = rm_CC(epsevTR,epsesTR,Ptr,Qtr,...
                                             epspvn,epspsn,kappa,mu0,alfa,P0,...
                                             epsev0,Pi0,A,pp,epstar,m,c0,mbar,Pin);
                                         
  end    
   
   if flag == 0 % ELASTIC STEP
  
        
    %Update variable
    epsev = epsevTR;
    epses = epsesTR;

    Pi= Pin;
    P = Ptr;
    Q = Qtr;
    epspv = epspvn;
    epsps = epspsn;

    % Compute matrix De 

    De = DE_HY(epsev,epses,P,kappa,mu0,alfa,P0,epsev0);

    Dep = De;
        
   end
       
  
OMEGA = -(epsev-epsev0)/kappa;
fac = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));  % Factor for computing alfa tensor in case of NaN

end

