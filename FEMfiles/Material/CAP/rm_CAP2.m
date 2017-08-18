function [P,Q,epsev,epses,epspv,epsps,Pi,Dep,flag] = ...
rm_CAP2(epsevTR,epsesTR,Ptr,Qtr,epspvn,epspsn,K,G,m,mbar,c0,A,pp,Pi0,epstar,Pin)

% Parameter
imax = 20;
toll = 10e-6;
NORMErec= zeros(imax);

% --- CAP 22 Return Mapping
RM = 1;   
  
% Inizialize variables
x = zeros(3,1);
x(1,1) = epsevTR;   % epsev
x(2,1) = epsesTR;   % epses
x(3,1) = 0;         % dgamma 

Pi = Pin;
P = Ptr;
Q = Qtr;

for iter = 1:imax

   % evaluate residual                
    B2 = m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2;
    F = B2*(P-Pi)^2+(A*Q)^2-B2*A^2;

    dPF = 2*B2*(P-Pi);
    dQF = 2*(A^2)*Q;

    r = [x(1) - epsevTR+x(3)*dPF
         x(2) - epsesTR+x(3)*dQF
               F                ];

       if iter == 1
       r0 = norm(r);
       end

    NORMErec(iter) = norm(r)/norm(r0);
    
    % check for convergence
    if norm(r) < toll*r0

        break
    else

        % evaluate tangent matrix
    Atang = Atang_CAP2(x(1),x(3),P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR);

        % solve for displacement increment
        dx = - (Atang\r);

        x = x + dx;                 
    end  

   [P,Q]=PQ_EL(x(1,1),x(2,1),K,G);

   epspv = epsevTR-x(1)+epspvn;
   epsps = epsesTR-x(2)+epspsn;

   Pi = Pi0*(epstar/(epstar-epspv))^pp;    
 
end

if iter == imax
    fprintf('No convergence of CAP2-RM')
%    NORMErec
end

% Compute ptilde
b = sqrt(m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2); 
ptilde = Pi-(mbar*A*b)/(sqrt((b^4/A^2)+mbar^2*b^2)); 


% --- CAP 21 Return Mapping
if P > ptilde  

RM = 2;        

% Inizialize variables
x = zeros(3,1);
x(1,1) = epsevTR;   % epsev
x(2,1) = epsesTR;   % epses
x(3,1) = 0;         % dgamma 

Pi = Pin;
P = Ptr;
Q = Qtr;   
    
    for iter = 1:imax 
        
    % evaluate residual                
    B2 = m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2;
    F = B2*(P-Pi)^2+(A*Q)^2-B2*A^2;
        
    dPG = -mbar;
    dQG = 1;

    r = [x(1) - epsevTR+x(3)*dPG
         x(2) - epsesTR+x(3)*dQG
               F                ];

       if iter == 1
       r0 = norm(r);
       end

    NORMErec(iter) = norm(r)/r0;

    % check for convergence
    if norm(r) < toll*r0

        break
    else

        % evaluate tangent matrix
    Atang = Atang_CAP21( x(1),P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR);

        % solve for displacement increment
        dx = - (Atang\r);

        x = x + dx;  

    end

    % Update

        [P,Q]=PQ_EL(x(1,1),x(2,1),K,G);

        epspv = epsevTR-x(1)+epspvn;
        epsps = epsesTR-x(2)+epspsn;


        Pi = Pi0*(epstar/(epstar-epspv))^pp;              

    end
end 

if iter == imax
    fprintf('No convergence of CAP21-RM')
%    NORMErec
end


%Update variable
epsev = x(1,1);
epses = x(2,1);
dgamma = x(3,1);

% Compute matrix Dep  
if RM == 1;
Dep = DEP_CAP2(epsev,dgamma,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR);
flag = 22;
else
Dep = DEP_CAP21(epsev,P,Q,Pi,K,G,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR);
flag = 21;
end

end

