function [P,Q,Dep,epsev,epses,epspv,epsps,Pi,flag,fac,NORMErec] = tens_CAP(Ge,epsevTR,epsesTR,epspvn,epspsn,Pin)

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
F_toll = 10e-8; % Attenzione: la tolleranza deve essere coerente 
                % con la grandezza di PQ del problema
NORMErec = zeros(15,1);
% Set isotropic elasto-plastic parameter

% 3	4	5  6	7	 8	  9    10	  11	
% E ni  m  Pi0  A    c0   pp  epstar  mbar

E = Ge(3);  
ni = Ge(4);
m = Ge(5);
Pi0 = Ge(6);
A = Ge(7);
c0 = Ge(8);
pp = Ge(9);
epstar = Ge(10);
mbar = Ge(11);

G = E/(2*(1+ni));
K = E/(3*(1-2*ni));

% Compute TRial stress    

[Ptr,Qtr] = PQ_EL(epsevTR,epsesTR,K,G);

flag = 0;

% Compute trial yield surface 

    % CAP2 yield surface
    B2 = m^2*Pin^2-m^2*A^2+2*m*c0*Pin+c0^2;
    F2tr = B2*(Ptr-Pin)^2+(A^2)*(Qtr^2)-(A^2)*B2;

    % CAP1 yield surface
    F1tr = Qtr-m*Ptr-c0;

    % Intersection CAP21
    b = sqrt(B2);
    pstar = Pin-(m*A*b)/(sqrt((b^4/A^2)+m^2*b^2));
           
% Check for plasticity

   if F2tr > F_toll 
  
      if F1tr > F_toll 
         flag = 1;    
         surf = 1; 
        
         % ------- Linear CAP_1 rm ------- %
       
         [P,Q,epsev,epses,epspv,epsps,Pi,Dep] = rm_CAP1(epsevTR,epsesTR,Ptr,Qtr,...
                                                       epspvn,epspsn,K,G,m,c0,mbar,Pin);
            
         % Check for apex RM
         if Q < 0
         fprintf('Apex region RM')
         end

         if P < pstar
           surf = 2; % Goto Cap2 surface
         end
                
      else
         if Ptr < pstar
            flag = 1;
            surf = 2;   
         end         
      end

   end
       
  if flag ==1 && surf == 2
         
[P,Q,epsev,epses,epspv,epsps,Pi,Dep,flag] = ...
rm_CAP2(epsevTR,epsesTR,Ptr,Qtr,epspvn,epspsn,K,G,m,mbar,c0,A,pp,Pi0,epstar,Pin);
        
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

    De = DE_EL(K,G);

    Dep = De;
        
end
       

fac = G;  % Factor for computing alfa tensor in case of NaN

end

